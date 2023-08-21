/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Edited by Thomas Campbell (Oxford)
------------------------------------------------------------------------- */

#include "pair_coul_cut_SPH.h"
#include <mpi.h>
#include <cmath>
#include <cstring>
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include "utils.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairCoulCutSPH::PairCoulCutSPH(LAMMPS *lmp) : Pair(lmp) {
  nmax = 0;

  centroidstressflag = 1;
  manybody_flag = 1;

  theta_coul = NULL;
  chi_coul_ei = NULL;

  comm_forward = 2;
  comm_reverse = 2;

}

/* ---------------------------------------------------------------------- */

PairCoulCutSPH::~PairCoulCutSPH()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(scale);

    memory->destroy(theta_coul);
    memory->destroy(chi_coul_ei);
  }
}

/* ---------------------------------------------------------------------- */

void PairCoulCutSPH::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,ecoul,fpair;
  double rsq,rsqrt,r2inv,rinv,forcecoul,factor_coul;
  int *ilist,*jlist,*numneigh,**firstneigh;

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

  double rho_i, omega_i, fact_i;

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
            theta_coul[i] += factor_coul*h2_i*scale[itype][jtype]*fact_i*q[j]*exp(-rsq/(2*eff_width*eff_width))/(eff_width*eff_width*eff_width);
          
            if (newton_pair || j < nlocal) {
              theta_coul[j] += factor_coul*scale[itype][jtype]*theta_const*((qtmp*q[j])/(rho_SPH[j]*omega_SPH[j]))*(width_SPH[j]*width_SPH[j]/(eff_width*eff_width*eff_width))*exp(-rsq/(2*eff_width*eff_width));
            }
          }

          if (jtype == ion_species){
            // ion neighbour
            chi_coul_ei[i] += factor_coul*scale[itype][jtype]*fact_i*q[j]*exp(-rsq/(2*h2_i))/(h_i);
          }
        }

        if (itype == ion_species) {
          // ion target
          if (jtype != ion_species){
            // electron neighbour
            if (newton_pair || j < nlocal) {
              chi_coul_ei[j] += factor_coul*scale[itype][jtype]*theta_const*((qtmp*q[j])/(rho_SPH[j]*omega_SPH[j]))*(1/width_SPH[j])*exp(-rsq/(2*width_SPH[j]*width_SPH[j]));
            }
          }
        }
      }
    }
  }

  if (newton_pair) comm->reverse_comm_pair(this);
  comm->forward_comm_pair(this);

  // loop over neighbors of my atoms

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
      // electron target

      h_i = width_SPH[i];
      h2_i = h_i*h_i;
      hm2_i = 1/h2_i;

      gauss_pre_i = pi_fact*(1/(h_i*h_i*h_i));
      // ele-ele and ion-ele SPH dynamic terms
      f[i][0] += -(theta_coul[i]+chi_coul_ei[i])*dx_rho_SPH[i];
      f[i][1] += -(theta_coul[i]+chi_coul_ei[i])*dy_rho_SPH[i];
      f[i][2] += -(theta_coul[i]+chi_coul_ei[i])*dz_rho_SPH[i];
   
    }
    
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_coul = special_coul[sbmask(j)];
    
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      rsqrt = pow(rsq,0.5);
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {

        if (jtype != ion_species){

          // electron neighbour

          h_j = width_SPH[j];
          hm2_j = 1/(h_j*h_j);

          if (itype == ion_species){

            // ion target

            force_fact = factor_coul*qqrd2e*scale[itype][jtype]*(qtmp*q[j])*(erf(rsqrt/(sqrt2*h_j))/(rsqrt*rsqrt*rsqrt) - (sqrt2/sqrt_pi)*(exp(-rsq*hm2_j/2)/(h_j*rsq)));

            f[i][0] += delx*force_fact;
            f[i][1] += dely*force_fact;
            f[i][2] += delz*force_fact;

            if (newton_pair || j < nlocal) {
              f[j][0] -= delx*force_fact;
              f[j][1] -= dely*force_fact;
              f[j][2] -= delz*force_fact;

            }
            if (eflag)ecoul = factor_coul * qqrd2e * scale[itype][jtype] * qtmp*q[j]*erf(rsqrt/(sqrt2*h_j))/rsqrt;
            if (evflag) ev_tally(i,j,nlocal,newton_pair,
                                 0.0,ecoul,force_fact,delx,dely,delz);
          }

          if (itype != ion_species){

            // electron target

            eff_width = pow((h2_i + width_SPH[j]*width_SPH[j]),0.5);

            force_fact = factor_coul*qqrd2e*scale[itype][jtype]*(qtmp*q[j])*(erf(rsqrt/(sqrt2*eff_width))/(rsqrt*rsqrt*rsqrt) - (sqrt2/sqrt_pi)*(exp(-rsq/(2*eff_width*eff_width))/(eff_width*rsq)));

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

            
            if (newton_pair || j < nlocal) {

              f[j][0] -= delx*force_fact;
              f[j][1] -= dely*force_fact;
              f[j][2] -= delz*force_fact;

              jmass = mass[jtype];
              m_gauss_ij = jmass*gauss_pre_i*exp(-(rsq)*hm2_i/2);

              ij_fact = hm2_i*m_gauss_ij*(theta_coul[i]+chi_coul_ei[i]);
              
              f[j][0] -= (delx)*ij_fact;
              f[j][1] -= (dely)*ij_fact;
              f[j][2] -= (delz)*ij_fact;

            }

            if (eflag) ecoul = factor_coul * qqrd2e*scale[itype][jtype]*qtmp*q[j]*erf(rsqrt/(sqrt2*eff_width))/rsqrt;
            // dynamic coulomb-SPH force expression is not pairwise symmetric
            // use of ev_tally not accurate for pressure evaluation - edit in future.
            full_factor = force_fact + 0.5*(ij_fact+ji_fact);
            if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,0.0,ecoul,
                        full_factor*delx - 0.5*(theta_coul[i]+chi_coul_ei[i])*dx_rho_SPH[i] + 0.5*(theta_coul[j]+chi_coul_ei[j])*dx_rho_SPH[j],full_factor*dely - 0.5*(theta_coul[i]+chi_coul_ei[i])*dy_rho_SPH[i] + 0.5*(theta_coul[j]+chi_coul_ei[j])*dy_rho_SPH[j],full_factor*delz - 0.5*(theta_coul[i]+chi_coul_ei[i])*dz_rho_SPH[i] + 0.5*(theta_coul[j]+chi_coul_ei[j])*dz_rho_SPH[j],delx,dely,delz);
          }
        }

        if (jtype == ion_species){

          // fprintf(screen,"\nin ion neighbour condition ");


          // ion neighbour

          if (itype != ion_species){

            // electron target

            force_fact = factor_coul*qqrd2e*scale[itype][jtype]*(qtmp*q[j])*(erf(rsqrt/(sqrt2*h_i))/(rsqrt*rsqrt*rsqrt) - (sqrt2/sqrt_pi)*(exp(-rsq*hm2_i/2)/(h_i*rsq)));

            f[i][0] += delx*force_fact;
            f[i][1] += dely*force_fact;
            f[i][2] += delz*force_fact;

            if (newton_pair || j < nlocal) {
              f[j][0] -= delx*force_fact;
              f[j][1] -= dely*force_fact;
              f[j][2] -= delz*force_fact;

            }
            if (eflag) ecoul = factor_coul * qqrd2e * scale[itype][jtype] * qtmp*q[j]*erf(rsqrt/(sqrt2*h_i))/rsqrt;
            if (evflag) ev_tally(i,j,nlocal,newton_pair,
                                 0.0,ecoul,force_fact,delx,dely,delz);
          }

          if (itype == ion_species){

            // ion target

            r2inv = 1.0/rsq;
            rinv = sqrt(r2inv);
            forcecoul = qqrd2e * scale[itype][jtype] * qtmp*q[j]*rinv;
            fpair = factor_coul * forcecoul * r2inv;

            f[i][0] += delx*fpair;
            f[i][1] += dely*fpair;
            f[i][2] += delz*fpair;

            if (newton_pair || j < nlocal) {
              f[j][0] -= delx*fpair;
              f[j][1] -= dely*fpair;
              f[j][2] -= delz*fpair;

            }
            if (eflag) ecoul = factor_coul * qqrd2e * scale[itype][jtype] * qtmp*q[j]*rinv;
          
            if (evflag) ev_tally(i,j,nlocal,newton_pair,
                                0.0,ecoul,fpair,delx,dely,delz);
          }
        }
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();

}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairCoulCutSPH::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(scale,n+1,n+1,"pair:scale");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairCoulCutSPH::settings(int narg, char **arg)
{
  if (narg != 5) error->all(FLERR,"Illegal pair_style command");

  cut_global = force->numeric(FLERR,arg[0]);
  ke_in = force->numeric(FLERR,arg[1]);
  ion_species = force->numeric(FLERR,arg[2]);
  pseudo_key = force->numeric(FLERR,arg[3]);
  ion_charge = force->numeric(FLERR,arg[4]);


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
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairCoulCutSPH::coeff(int narg, char **arg)
{
  if (narg < 2 || narg > 3)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  double cut_one = cut_global;
  if (narg == 3) cut_one = force->numeric(FLERR,arg[2]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      cut[i][j] = cut_one;
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

void PairCoulCutSPH::init_style()
{
  if (!atom->q_flag)
    error->all(FLERR,"Pair style coul/cut requires atom attribute q");

  neighbor->request(this,instance_me);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairCoulCutSPH::init_one(int i, int j)
{
  if (setflag[i][j] == 0)
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);

  scale[j][i] = scale[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairCoulCutSPH::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) fwrite(&cut[i][j],sizeof(double),1,fp);
    }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairCoulCutSPH::read_restart(FILE *fp)
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
        if (me == 0) utils::sfread(FLERR,&cut[i][j],sizeof(double),1,fp,NULL,error);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairCoulCutSPH::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairCoulCutSPH::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR,&cut_global,sizeof(double),1,fp,NULL,error);
    utils::sfread(FLERR,&offset_flag,sizeof(int),1,fp,NULL,error);
    utils::sfread(FLERR,&mix_flag,sizeof(int),1,fp,NULL,error);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ---------------------------------------------------------------------- */

void *PairCoulCutSPH::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str,"cut_coul") == 0) return (void *) &cut;
  if (strcmp(str,"scale") == 0) return (void *) scale;
  return NULL;
}

/* ---------------------------------------------------------------------- */


int PairCoulCutSPH::pack_forward_comm(int n, int *list, double *buf,
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

void PairCoulCutSPH::unpack_forward_comm(int n, int first, double *buf)
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

int PairCoulCutSPH::pack_reverse_comm(int n, int first, double *buf)
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

void PairCoulCutSPH::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    theta_coul[j] += buf[m++];
    chi_coul_ei[j] += buf[m++];
  }
  
}
