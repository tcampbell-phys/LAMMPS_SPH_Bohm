/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Paul Crozier (SNL)
   Edited by Thomas Campbell (Oxford)
   - for use with dynamic per-particle SPH kernel widths
   - 'PH' suffix refers to 'post-hoc' removal of erf contributions to Coulomb interactions
     in the real space calculations.
------------------------------------------------------------------------- */

#include "pair_coul_long_SPH_ph.h"
#include "domain.h"
#include <mpi.h>
#include <cmath>
#include <cstring>
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "kspace.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include "utils.h"

using namespace LAMMPS_NS;

#define EWALD_F   1.12837917
#define EWALD_P   0.3275911
#define A1        0.254829592
#define A2       -0.284496736
#define A3        1.421413741
#define A4       -1.453152027
#define A5        1.061405429

/* ---------------------------------------------------------------------- */

PairCoulLongSPHph::PairCoulLongSPHph(LAMMPS *lmp) : Pair(lmp)
{
  ewaldflag = pppmflag = 1;
  ftable = NULL;
  qdist = 0.0;
  cut_respa = NULL;
  manybody_flag = 1;

  nmax = 0;

  numforce = 0;

  theta_coul = NULL;
  chi_coul_ei = NULL;

  comm_forward = 2;
  comm_reverse = 2;
}

/* ---------------------------------------------------------------------- */

PairCoulLongSPHph::~PairCoulLongSPHph()
{
  if (copymode) return;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(theta_coul);
    memory->destroy(chi_coul_ei);

    memory->destroy(scale);
  }
  // still require this condition without use of tables in compute?
  if (ftable) free_tables();
}

/* ---------------------------------------------------------------------- */

void PairCoulLongSPHph::compute(int eflag, int vflag)
{
  int i,j,ii,jj,k,inum,jnum,itable,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,ecoul,fpair,ecoul_erfc,ecoul_erf,fpair_erf;
  double fraction,table,energy_table;
  double r,r2inv,forcecoul,factor_coul,coul_prefact;
  double grij,expm2,prefactor,t,erfc;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double rsq;

  ecoul = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  double *special_coul = force->special_coul;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;
  double *rho_SPH = atom->rho_SPH;
  double *width_SPH = atom->width_SPH;
  double *omega_SPH = atom->omega_SPH;
  double *dx_rho_SPH = atom->dx_rho_SPH;
  double *dy_rho_SPH = atom->dy_rho_SPH;
  double *dz_rho_SPH = atom->dz_rho_SPH;

  double rho_i, omega_i, fact_i, fact_chi, mu_k, mu_k_min3_2, beta_k, prefact;

  double eff_width;
  double force_fact;

  double gauss_pre_j,gauss_pre_i;

  double h_i,h2_i,hm2_i,hm4_i;
  double h_j,h2_j,hm2_j,hm4_j;
  double imass,jmass,ijmass;
  double m_gauss_ij,m_gauss_ji;
  double ij_fact,ji_fact;
  double full_factor;
  double *mass = atom->mass;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  if (atom->nmax > nmax) {
    // delete and create new memory arrays for any per-particle variables that need communicating.
    memory->destroy(theta_coul);
    memory->destroy(chi_coul_ei);
    nmax = atom->nmax;
    memory->create(theta_coul,nmax,"pair:theta_coul");
    memory->create(chi_coul_ei,nmax,"pair:chi_coul_ei");
  }
  // zero out per-atom arrays

  if (newton_pair) {
    for (i = 0; i < nall; i++){
      theta_coul[i] = 0.0;
      chi_coul_ei[i] = 0.0;
    }
  } 
  else{
    for (i = 0; i < nlocal; i++){
      theta_coul[i] = 0.0;
      chi_coul_ei[i] = 0.0;
    }
  }

  double theta_const;

  theta_const = qqrd2e*sqrt2/(sqrt_pi*3);

  double pi_fact = 1/(sqrt2*sqrt_pi*sqrt2*sqrt_pi*sqrt2*sqrt_pi);

  // compute communicated properties

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    imass = mass[itype];
    if (itype != ion_species) {
      h_i = width_SPH[i];
      h2_i = h_i*h_i;

      fact_i = (theta_const*qtmp)/(rho_SPH[i]*omega_SPH[i]);
      fact_chi = h_i/(3*rho_SPH[i]*omega_SPH[i]);
    }

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_coul = special_coul[sbmask(j)];
      
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {

        if (itype != ion_species) {
          // electron target
          if (jtype != ion_species){
            // electron neighbour
            eff_width = pow((h2_i + width_SPH[j]*width_SPH[j]),0.5);
            theta_coul[i] += h2_i*scale[itype][jtype]*fact_i*q[j]*exp(-rsq/(2*eff_width*eff_width))/(eff_width*eff_width*eff_width);
          
            if (newton_pair || j < nlocal) {
              theta_coul[j] += scale[itype][jtype]*theta_const*((qtmp*q[j])/(rho_SPH[j]*omega_SPH[j]))*(width_SPH[j]*width_SPH[j]/(eff_width*eff_width*eff_width))*exp(-rsq/(2*eff_width*eff_width));
            }
          }

          if (jtype == ion_species){
            // ion neighbour
            prefact = scale[itype][jtype]*fact_chi;
             for (k = 0; k < N_coeff; k++){
              mu_k = 2*h2_i*a_coeff[k] + 1;
              beta_k = a_coeff[k]/mu_k;
              chi_coul_ei[i] += prefact*c_coeff[k]*2*h_i*a_coeff[k]*exp(-beta_k*rsq)*pow(mu_k,-2.5)*(2*rsq*beta_k - 3);
            }
          }
        }

        if (itype == ion_species) {
          // ion target
          if (jtype != ion_species){
            // electron neighbour
            if (newton_pair || j < nlocal) {
              prefact = scale[itype][jtype]*(width_SPH[j]/(3*rho_SPH[j]*omega_SPH[j]));
              for (k = 0; k < N_coeff; k++){
                mu_k = 2*width_SPH[j]*width_SPH[j]*a_coeff[k] + 1;
                beta_k = a_coeff[k]/mu_k;
                chi_coul_ei[j] += prefact*c_coeff[k]*2*width_SPH[j]*a_coeff[k]*exp(-beta_k*rsq)*pow(mu_k,-2.5)*(2*rsq*beta_k - 3);
              }
            }
          }
        }
      }
    }
  }

  if (newton_pair) comm->reverse_comm_pair(this);
  comm->forward_comm_pair(this);

  // loop over neighbors of my atoms

  double e_fx_dyn, e_fy_dyn, e_fz_dyn;
  double ie_fx, ie_fy, ie_fz;
  double ei_fx, ei_fy, ei_fz;
  double ee_ij_fx, ee_ij_fy, ee_ij_fz;
  double ee_ji_fx, ee_ji_fy, ee_ji_fz;

  double force_limit = 1000;
  double dyn_fact_limit = 4000000;


  for (ii = 0; ii < inum; ii++) {

    i = ilist[ii];
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    imass = mass[itype];

    if (itype != ion_species) {
      // fprintf(screen,"\n\ne target...");
      // electron target

      h_i = width_SPH[i];
      h2_i = h_i*h_i;
      hm2_i = 1/h2_i;

      gauss_pre_i = pi_fact*(1/(h_i*h_i*h_i));
      // fprintf(screen,"\ne fx self = %16.16f",-(theta_coul[i]+chi_coul_ei[i])*dx_rho_SPH[i]);
      // fprintf(screen,"\ne fy self = %16.16f",-(theta_coul[i]+chi_coul_ei[i])*dy_rho_SPH[i]);
      // fprintf(screen,"\ne fz self = %16.16f",-(theta_coul[i]+chi_coul_ei[i])*dz_rho_SPH[i]);

      e_fx_dyn = -(theta_coul[i]+chi_coul_ei[i])*dx_rho_SPH[i];
      e_fy_dyn = -(theta_coul[i]+chi_coul_ei[i])*dy_rho_SPH[i];
      e_fz_dyn = -(theta_coul[i]+chi_coul_ei[i])*dz_rho_SPH[i];

      // if (e_fx_dyn > force_limit || e_fy_dyn > force_limit || e_fz_dyn > force_limit){
      //   fprintf(screen,"\n\nWARNING: large force in e dyn...");
      //   fprintf(screen,"\ne x = %16.16f",xtmp);
      //   fprintf(screen,"\ne y = %16.16f",ytmp);
      //   fprintf(screen,"\ne z = %16.16f",ztmp);
      //   fprintf(screen,"\ne_fx_dyn = %16.16f",e_fx_dyn);
      //   fprintf(screen,"\ne_fy_dyn = %16.16f",e_fy_dyn);
      //   fprintf(screen,"\ne_fz_dyn = %16.16f",e_fz_dyn);
      // }

      // ele-ele and ion-ele SPH dynamic terms
      f[i][0] += -(theta_coul[i]+chi_coul_ei[i])*dx_rho_SPH[i];
      f[i][1] += -(theta_coul[i]+chi_coul_ei[i])*dy_rho_SPH[i];
      f[i][2] += -(theta_coul[i]+chi_coul_ei[i])*dz_rho_SPH[i];

      // if (theta_coul[i] > dyn_fact_limit || chi_coul_ei[i] > dyn_fact_limit){
      //   fprintf(screen,"\n\nWARNING: large dynamic electron factor in e dyn...");
      //   fprintf(screen,"\ne x = %16.16f",xtmp);
      //   fprintf(screen,"\ne y = %16.16f",ytmp);
      //   fprintf(screen,"\ne z = %16.16f",ztmp);
      //   fprintf(screen,"\ntheta_coul[i] = %16.16f",theta_coul[i]);
      //   fprintf(screen,"\nchi_coul_ei[i] = %16.16f",chi_coul_ei[i]);
      // }
    
    }
    else{
      // fprintf(screen,"\n\ni target...");
    }


    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      factor_coul = special_coul[sbmask(j)];

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cut_coulsq) {
        // fprintf(screen,"\nrsq = %16.16f",rsq);

        coul_prefact = qqrd2e*scale[itype][jtype]*qtmp*q[j];

        r2inv = 1.0/rsq;
        r = sqrt(rsq);
        // // fprintf(screen,"\n\ng = %16.16f",g_ewald);
        // // fprintf(screen,"\nr = %16.16f",r);
        // // fprintf(screen,"\ncharge_prod = %16.16f",qtmp*q[j]);
        prefactor = coul_prefact/r;
        //use standard coul/long force loop to extract erfc and hence erf terms
        if (!ncoultablebits || rsq <= tabinnersq) {
          grij = g_ewald * r;
          expm2 = exp(-grij*grij);
          t = 1.0 / (1.0 + EWALD_P*grij);
          erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;
          ecoul_erfc = prefactor*erfc;
          forcecoul = prefactor * (erfc + EWALD_F*grij*expm2);
        } else {
          union_int_float_t rsq_lookup;
          rsq_lookup.f = rsq;
          itable = rsq_lookup.i & ncoulmask;
          itable >>= ncoulshiftbits;
          fraction = (rsq_lookup.f - rtable[itable]) * drtable[itable];
          table = ftable[itable] + fraction*dftable[itable];
          forcecoul = scale[itype][jtype] * qtmp*q[j] * table;
          energy_table = etable[itable] + fraction*detable[itable];
          ecoul_erfc = scale[itype][jtype] * qtmp*q[j] * energy_table;
        }

        fpair = forcecoul * r2inv;

        fpair_erf = -fpair + prefactor*r2inv;
        ecoul_erf = -ecoul_erfc + prefactor;

        // // fprintf(screen,"\nfpair erfc = %16.16f",fpair);
        // // fprintf(screen,"\nfpair erf = %16.16f",fpair_erf);
        // // fprintf(screen,"\nfpair = %16.16f",prefactor*r2inv);

        // // fprintf(screen,"\necoul erfc = %16.16f",ecoul_erfc);
        // // fprintf(screen,"\necoul erf = %16.16f",ecoul_erf);
        // // fprintf(screen,"\necoul = %16.16f",prefactor);

        if (jtype != ion_species){

          // electron neighbour

          h_j = width_SPH[j];
          hm2_j = 1/(h_j*h_j);

          if (itype == ion_species){

            // ion target
            ecoul = 0.;
            force_fact = 0.;

            for (k = 0; k < N_coeff; k++){
              mu_k = 2*h_j*h_j*a_coeff[k] + 1;
              beta_k = a_coeff[k]/mu_k;
              mu_k_min3_2 = pow(mu_k,-1.5);
              ecoul += c_coeff[k] * mu_k_min3_2 * exp(-beta_k*rsq);
              force_fact += 2 * c_coeff[k] * beta_k * mu_k_min3_2 * exp(-beta_k*rsq);
            }

            ecoul *=  qqrd2e * scale[itype][jtype];

            // remove erf coul

            ecoul -= ecoul_erf;

            force_fact *=  qqrd2e * scale[itype][jtype];

            // ie_fx = delx*force_fact;
            // ie_fy = dely*force_fact;
            // ie_fz = delz*force_fact;

            // if (ie_fx > force_limit || ie_fy > force_limit || ie_fz > force_limit){
            //   fprintf(screen,"\n\nWARNING: large force in i-e interaction...");
            //   fprintf(screen,"\ni x = %16.16f",xtmp);
            //   fprintf(screen,"\ni y = %16.16f",ytmp);
            //   fprintf(screen,"\ni z = %16.16f",ztmp);
            //   fprintf(screen,"\ne x = %16.16f",x[j][0]);
            //   fprintf(screen,"\ne y = %16.16f",x[j][1]);
            //   fprintf(screen,"\ne z = %16.16f",x[j][2]);
            //   fprintf(screen,"\nie_fx = %16.16f",ie_fx);
            //   fprintf(screen,"\nie_fy = %16.16f",ie_fy);
            //   fprintf(screen,"\nie_fz = %16.16f",ie_fz);
            // }
            
            f[i][0] += delx*force_fact;
            f[i][1] += dely*force_fact;
            f[i][2] += delz*force_fact;

            // fprintf(screen,"\ni-e coul = %16.16f",ecoul);
            // fprintf(screen,"\ni-e fact = %16.16f",force_fact);

            // remove erf coul forces

            f[i][0] -= (delx)*fpair_erf;
            f[i][1] -= (dely)*fpair_erf;
            f[i][2] -= (delz)*fpair_erf;

            // fprintf(screen,"\ni-e fx = %16.16f",delx*force_fact);
            // fprintf(screen,"\ni-e fy = %16.16f",dely*force_fact);
            // fprintf(screen,"\ni-e fz = %16.16f",delz*force_fact);


            if (newton_pair || j < nlocal) {
              f[j][0] -= delx*force_fact;
              f[j][1] -= dely*force_fact;
              f[j][2] -= delz*force_fact;

              // fprintf(screen,"\ne-i fx = %16.16f",-delx*force_fact);
              // fprintf(screen,"\ne-i fy = %16.16f",-dely*force_fact);
              // fprintf(screen,"\ne-i fz = %16.16f",-delz*force_fact);


              f[j][0] += delx*fpair_erf;
              f[j][1] += dely*fpair_erf;
              f[j][2] += delz*fpair_erf;
            }

            if (evflag) ev_tally(i,j,nlocal,newton_pair,
                                 0.0,ecoul,force_fact-fpair_erf,delx,dely,delz);
          }

          if (itype != ion_species){

            // electron target

            // ecoul = qqrd2e*scale[itype][jtype]*qtmp*q[j]*erf(r/(sqrt2*eff_width))/r;
            
            // remove erf coul

            ecoul -= ecoul_erf;

            eff_width = pow((h2_i + width_SPH[j]*width_SPH[j]),0.5);

            force_fact = coul_prefact*(erf(r/(sqrt2*eff_width))/(r*r*r) - (sqrt2/sqrt_pi)*(exp(-rsq/(2*eff_width*eff_width))/(eff_width*rsq)));
            
            f[i][0] += delx*force_fact;
            f[i][1] += dely*force_fact;
            f[i][2] += delz*force_fact;

            h_j = width_SPH[j];
            hm2_j = 1/(h_j*h_j);
            gauss_pre_j = pi_fact*(1/(h_j*h_j*h_j));
            m_gauss_ji = imass*gauss_pre_j*exp(-(rsq)*hm2_j/2);
            ji_fact = hm2_j*m_gauss_ji*(theta_coul[j]+chi_coul_ei[j]);

            // ele-ele and ion-ele SPH dynamic width terms

            f[i][0] += (delx)*ji_fact;
            f[i][1] += (dely)*ji_fact;
            f[i][2] += (delz)*ji_fact;

            // fprintf(screen,"\ne-e ji fx = %16.16f",delx*(force_fact+ji_fact));
            // fprintf(screen,"\ne-e ji fy = %16.16f",dely*(force_fact+ji_fact));
            // fprintf(screen,"\ne-e ji fz = %16.16f",delz*(force_fact+ji_fact));

            // ee_ji_fx = delx*(force_fact+ji_fact);
            // ee_ji_fy = dely*(force_fact+ji_fact);
            // ee_ji_fz = delz*(force_fact+ji_fact);

            // if (ee_ji_fx > force_limit || ee_ji_fy > force_limit || ee_ji_fz > force_limit){
            //   fprintf(screen,"\n\nWARNING: large force in e-e ji interaction...");
            //   fprintf(screen,"\ne i x = %16.16f",xtmp);
            //   fprintf(screen,"\ne i y = %16.16f",ytmp);
            //   fprintf(screen,"\ne i z = %16.16f",ztmp);
            //   fprintf(screen,"\ne j x = %16.16f",x[j][0]);
            //   fprintf(screen,"\ne j y = %16.16f",x[j][1]);
            //   fprintf(screen,"\ne j z = %16.16f",x[j][2]);
            //   fprintf(screen,"\ne-e rsq = %16.16f",rsq);
            //   fprintf(screen,"\nh_i = %16.16f",h_i);
            //   fprintf(screen,"\nh_j = %16.16f",h_j);
            //   fprintf(screen,"\nee_ji_fx = %16.16f",ee_ji_fx);
            //   fprintf(screen,"\nee_ji_fy = %16.16f",ee_ji_fy);
            //   fprintf(screen,"\nee_ji_fz = %16.16f",ee_ji_fz);
            //   fprintf(screen,"\nee force_fact = %16.16f",force_fact);
            //   fprintf(screen,"\nee ji_fact = %16.16f",ji_fact);
            //   fprintf(screen,"\nee m_gauss_ji = %16.16f",m_gauss_ji);
            //   fprintf(screen,"\nee theta_coul[j] = %16.16f",theta_coul[j]);
            //   fprintf(screen,"\nee chi_coul_ei[j] = %16.16f",chi_coul_ei[j]);
            // }

            // remove erf coul forces

            f[i][0] -= (delx)*fpair_erf;
            f[i][1] -= (dely)*fpair_erf;
            f[i][2] -= (delz)*fpair_erf;

            if (newton_pair || j < nlocal) {

              f[j][0] -= delx*force_fact;
              f[j][1] -= dely*force_fact;
              f[j][2] -= delz*force_fact;

              f[j][0] += delx*fpair_erf;
              f[j][1] += dely*fpair_erf;
              f[j][2] += delz*fpair_erf;

              jmass = mass[jtype];
              m_gauss_ij = jmass*gauss_pre_i*exp(-(rsq)*hm2_i/2);

              ij_fact = hm2_i*m_gauss_ij*(theta_coul[i]+chi_coul_ei[i]);

              // ee_ij_fx = delx*(force_fact+ij_fact);
              // ee_ij_fy = dely*(force_fact+ij_fact);
              // ee_ij_fz = delz*(force_fact+ij_fact);

              // if (ee_ij_fx > force_limit || ee_ij_fy > force_limit || ee_ij_fz > force_limit){
              //   fprintf(screen,"\n\nWARNING: large force in e-e ij (newton) interaction...");
              //   fprintf(screen,"\ne i x = %16.16f",xtmp);
              //   fprintf(screen,"\ne i y = %16.16f",ytmp);
              //   fprintf(screen,"\ne i z = %16.16f",ztmp);
              //   fprintf(screen,"\ne j x = %16.16f",x[j][0]);
              //   fprintf(screen,"\ne j y = %16.16f",x[j][1]);
              //   fprintf(screen,"\ne j z = %16.16f",x[j][2]);
              //   fprintf(screen,"\nh_i = %16.16f",h_i);
              //   fprintf(screen,"\nh_j = %16.16f",h_j);
              //   fprintf(screen,"\nee_ij_fx = %16.16f",ee_ij_fx);
              //   fprintf(screen,"\nee_ij_fy = %16.16f",ee_ij_fy);
              //   fprintf(screen,"\nee_ij_fz = %16.16f",ee_ij_fz);
              //   fprintf(screen,"\nee force_fact = %16.16f",force_fact);
              //   fprintf(screen,"\nee ij_fact = %16.16f",ij_fact);
              // }
              
              f[j][0] -= (delx)*ij_fact;
              f[j][1] -= (dely)*ij_fact;
              f[j][2] -= (delz)*ij_fact;

              // fprintf(screen,"\ne-e ij fx = %16.16f",-delx*(force_fact+ij_fact));
              // fprintf(screen,"\ne-e ij fy = %16.16f",-dely*(force_fact+ij_fact));
              // fprintf(screen,"\ne-e ij fz = %16.16f",-delz*(force_fact+ij_fact));


            }

            if (eflag) ecoul = factor_coul * qqrd2e*scale[itype][jtype]*qtmp*q[j]*erf(r/(sqrt2*eff_width))/r;
            

            // fprintf(screen,"\ne-e coul = %16.16f",ecoul);
            // fprintf(screen,"\ne-e fact = %16.16f",force_fact + 0.5*(ij_fact+ji_fact));

            // dynamic coulomb-SPH force expression is not pairwise symmetric
            // use of ev_tally not accurate for pressure evaluation - edit in future.
            
            full_factor = force_fact - fpair_erf + 0.5*(ij_fact+ji_fact);
            if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,0.0,ecoul,
                        full_factor*delx - 0.5*(theta_coul[i]+chi_coul_ei[i])*dx_rho_SPH[i] + 0.5*(theta_coul[j]+chi_coul_ei[j])*dx_rho_SPH[j],full_factor*dely - 0.5*(theta_coul[i]+chi_coul_ei[i])*dy_rho_SPH[i] + 0.5*(theta_coul[j]+chi_coul_ei[j])*dy_rho_SPH[j],full_factor*delz - 0.5*(theta_coul[i]+chi_coul_ei[i])*dz_rho_SPH[i] + 0.5*(theta_coul[j]+chi_coul_ei[j])*dz_rho_SPH[j],delx,dely,delz);
          }
        }

        if (jtype == ion_species){

          // ion neighbour

          if (itype != ion_species){

            // electron target
            ecoul = 0.;
            force_fact = 0.;

            for (k = 0; k < N_coeff; k++){
              mu_k = 2*h2_i*a_coeff[k] + 1;
              beta_k = a_coeff[k]/mu_k;
              mu_k_min3_2 = pow(mu_k,-1.5);
              ecoul += c_coeff[k] * mu_k_min3_2 * exp(-beta_k*rsq);
              force_fact += 2 * c_coeff[k] * beta_k * mu_k_min3_2 * exp(-beta_k*rsq);
            }
            ecoul *=  qqrd2e * scale[itype][jtype];
            force_fact *=  qqrd2e * scale[itype][jtype];

            // fprintf(screen,"\ne-i coul = %16.16f",ecoul);
            // fprintf(screen,"\ne-i fact = %16.16f",force_fact);

            // remove erf coul

            ecoul -= ecoul_erf;
            
            f[i][0] += delx*force_fact;
            f[i][1] += dely*force_fact;
            f[i][2] += delz*force_fact;

            ei_fx = delx*force_fact;
            ei_fy = dely*force_fact;
            ei_fz = delz*force_fact;

            // if (ie_fx > force_limit || ie_fy > force_limit || ie_fz > force_limit){
            //   fprintf(screen,"\n\nWARNING: large force in e-i interaction...");
            //   fprintf(screen,"\ne x = %16.16f",xtmp);
            //   fprintf(screen,"\ne y = %16.16f",ytmp);
            //   fprintf(screen,"\ne z = %16.16f",ztmp);
            //   fprintf(screen,"\ni x = %16.16f",x[j][0]);
            //   fprintf(screen,"\ni y = %16.16f",x[j][1]);
            //   fprintf(screen,"\ni z = %16.16f",x[j][2]);
            //   fprintf(screen,"\nie_fx = %16.16f",ie_fx);
            //   fprintf(screen,"\nie_fy = %16.16f",ie_fy);
            //   fprintf(screen,"\nie_fz = %16.16f",ie_fz);
            // }

            // fprintf(screen,"\ne-i fx = %16.16f",delx*force_fact);
            // fprintf(screen,"\ne-i fy = %16.16f",dely*force_fact);
            // fprintf(screen,"\ne-i fz = %16.16f",delz*force_fact);


            // remove erf coul forces

            f[i][0] -= (delx)*fpair_erf;
            f[i][1] -= (dely)*fpair_erf;
            f[i][2] -= (delz)*fpair_erf;

            if (newton_pair || j < nlocal) {
              f[j][0] -= delx*force_fact;
              f[j][1] -= dely*force_fact;
              f[j][2] -= delz*force_fact;

              // fprintf(screen,"\ni-e fx = %16.16f",-delx*force_fact);
              // fprintf(screen,"\ni-e fy = %16.16f",-dely*force_fact);
              // fprintf(screen,"\ni-e fz = %16.16f",-delz*force_fact);


              f[j][0] += delx*fpair_erf;
              f[j][1] += dely*fpair_erf;
              f[j][2] += delz*fpair_erf;
            }
            if (evflag) ev_tally(i,j,nlocal,newton_pair,
                                 0.0,ecoul,force_fact-fpair_erf,delx,dely,delz);
          }

          if (itype == ion_species){

            // ion target

            // fprintf(screen,"\ni-i forcecoul = %16.16f",forcecoul);
            // fprintf(screen,"\ni-i factor_coul = %16.16f",factor_coul);

            f[i][0] += delx*fpair;
            f[i][1] += dely*fpair;
            f[i][2] += delz*fpair;

            // fprintf(screen,"\ni-i fx = %16.16f",delx*fpair);
            // fprintf(screen,"\ni-i fy = %16.16f",dely*fpair);
            // fprintf(screen,"\ni-i fz = %16.16f",delz*fpair);


            if (newton_pair || j < nlocal) {
              f[j][0] -= delx*fpair;
              f[j][1] -= dely*fpair;
              f[j][2] -= delz*fpair;


              // fprintf(screen,"\ni-i fx = %16.16f",-delx*fpair);
              // fprintf(screen,"\ni-i fy = %16.16f",-dely*fpair);
              // fprintf(screen,"\ni-i fz = %16.16f",-delz*fpair);


            }

            // fprintf(screen,"\ni-i coul = %16.16f",factor_coul * qqrd2e * scale[itype][jtype] * qtmp*q[j]/r);
            // fprintf(screen,"\ni-i fact = %16.16f",fpair);

            if (evflag) ev_tally(i,j,nlocal,newton_pair,
                                0.0,ecoul_erfc,fpair,delx,dely,delz);
          }
        }
      }
    }
  }
  numforce = 1;
  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairCoulLongSPHph::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(scale,n+1,n+1,"pair:scale");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairCoulLongSPHph::settings(int narg, char **arg)
{
  if (narg != 6) error->all(FLERR,"Illegal pair_style command");

  cut_coul = force->numeric(FLERR,arg[0]);
  ke_in = force->numeric(FLERR,arg[1]);
  ion_species = force->numeric(FLERR,arg[2]);
  pseudo_key = force->numeric(FLERR,arg[3]);
  ion_charge = force->numeric(FLERR,arg[4]);
  ele_TFWHM = force->numeric(FLERR,arg[5]);

  // check input parameters match selected Pseudopotential parameters
  if (pseudo_key==al_lda_A_key){
    c_coeff = al_lda_A_c;
    a_coeff = al_lda_A_a;
    N_coeff = al_lda_A_Ncoeff;
    if (ion_charge != al_lda_A_ion_charge){
      fprintf(screen,"\n ion_charge = %d \n",ion_charge);
      fprintf(screen,"\n al_lda_A_ion_charge = %d \n",al_lda_A_ion_charge);
      error->all(FLERR,"PSEUDO_ERR requested pseudopotential parameters do not match input: ion_charge ");
    }
    if (cut_coul < ele_TFWHM){
      fprintf(screen,"\n cut_coul = %16.16f \n",cut_coul);
      fprintf(screen,"\n ele_TFWHM = %16.16f \n",ele_TFWHM);
      error->all(FLERR,"PSEUDO_ERR cutoff too short for use with SPH coulomb: cutoff must be larger than TFWHM");
    }
  }
    if (pseudo_key==al_lda_D_key){
    c_coeff = al_lda_D_c;
    a_coeff = al_lda_D_a;
    N_coeff = al_lda_D_Ncoeff;
    if (ion_charge != al_lda_D_ion_charge){
      fprintf(screen,"\n ion_charge = %d \n",ion_charge);
      fprintf(screen,"\n al_lda_D_ion_charge = %d \n",al_lda_D_ion_charge);
      error->all(FLERR,"PSEUDO_ERR requested pseudopotential parameters do not match input: ion_charge ");
    }
    if (cut_coul < ele_TFWHM){
      fprintf(screen,"\n cut_coul = %16.16f \n",cut_coul);
      fprintf(screen,"\n ele_TFWHM = %16.16f \n",ele_TFWHM);
      error->all(FLERR,"PSEUDO_ERR cutoff too short for use with SPH coulomb: cutoff must be larger than TFWHM");
    }
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairCoulLongSPHph::coeff(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      scale[i][j] = 1.0;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairCoulLongSPHph::init_style()
{
  if (!atom->q_flag)
    error->all(FLERR,"Pair style lj/cut/coul/long requires atom attribute q");

  neighbor->request(this,instance_me);

  cut_coulsq = cut_coul * cut_coul;

  // ensure use of KSpace long-range solver, set g_ewald

 if (force->kspace == NULL)
    error->all(FLERR,"Pair style requires a KSpace style");
  g_ewald = force->kspace->g_ewald;

  // setup force tables

  if (ncoultablebits) init_tables(cut_coul,NULL);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairCoulLongSPHph::init_one(int i, int j)
{
  scale[j][i] = scale[i][j];
  return cut_coul+2.0*qdist;
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairCoulLongSPHph::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j])
        fwrite(&scale[i][j],sizeof(double),1,fp);
    }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairCoulLongSPHph::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR,&setflag[i][j],sizeof(int),1,fp,NULL,error);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) utils::sfread(FLERR,&scale[i][j],sizeof(double),1,fp,NULL,error);
        MPI_Bcast(&scale[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairCoulLongSPHph::write_restart_settings(FILE *fp)
{
  fwrite(&cut_coul,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&ncoultablebits,sizeof(int),1,fp);
  fwrite(&tabinner,sizeof(double),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairCoulLongSPHph::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR,&cut_coul,sizeof(double),1,fp,NULL,error);
    utils::sfread(FLERR,&offset_flag,sizeof(int),1,fp,NULL,error);
    utils::sfread(FLERR,&mix_flag,sizeof(int),1,fp,NULL,error);
    utils::sfread(FLERR,&ncoultablebits,sizeof(int),1,fp,NULL,error);
    utils::sfread(FLERR,&tabinner,sizeof(double),1,fp,NULL,error);
  }
  MPI_Bcast(&cut_coul,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
  MPI_Bcast(&ncoultablebits,1,MPI_INT,0,world);
  MPI_Bcast(&tabinner,1,MPI_DOUBLE,0,world);
}

/* ---------------------------------------------------------------------- */

// double PairCoulLongSPHph::single(int i, int j, int /*itype*/, int /*jtype*/,
//                             double rsq,
//                             double factor_coul, double /*factor_lj*/,
//                             double &fforce)
// {
//   double r2inv,r,grij,expm2,t,erfc,prefactor;
//   double fraction,table,forcecoul,phicoul;
//   int itable;

//   if (numforce==0){
//     error->all(FLERR,"coul/long/SPH many-boby data required for this calculation is missing");
//   }

//   if (numforce==1){

//     r2inv = 1.0/rsq;
//     r = sqrt(rsq);
//     grij = g_ewald * r;
//     expm2 = exp(-grij*grij);
//     t = 1.0 / (1.0 + EWALD_P*grij);
//     erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;
//     prefactor = force->qqrd2e * atom->q[i]*atom->q[j]/r;
//     forcecoul = prefactor * (erfc + EWALD_F*grij*expm2);
//     if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*prefactor;
    
//     fforce = forcecoul * r2inv;

//     phicoul = prefactor*erfc;
    
//     if (factor_coul < 1.0) phicoul -= (1.0-factor_coul)*prefactor;

//     return phicoul;
//   }
// }

/* ---------------------------------------------------------------------- */

void *PairCoulLongSPHph::extract(const char *str, int &dim)
{
  if (strcmp(str,"cut_coul") == 0) {
    dim = 0;
    return (void *) &cut_coul;
  }
  if (strcmp(str,"scale") == 0) {
    dim = 2;
    return (void *) scale;
  }
  return NULL;
}

/* ---------------------------------------------------------------------- */


int PairCoulLongSPHph::pack_forward_comm(int n, int *list, double *buf,
                               int /*pbc_flag*/, int * /*pbc*/)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = theta_coul[j];
    buf[m++] = chi_coul_ei[j];
  }
  return m;
  
}

/* ---------------------------------------------------------------------- */

void PairCoulLongSPHph::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++){
    theta_coul[i] = buf[m++];
    chi_coul_ei[i] = buf[m++];
  }
  
}
/* ---------------------------------------------------------------------- */

int PairCoulLongSPHph::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++){
    buf[m++] = theta_coul[i];
    buf[m++] = chi_coul_ei[i];
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void PairCoulLongSPHph::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    theta_coul[j] += buf[m++];
    chi_coul_ei[j] += buf[m++];
  }
  
}
