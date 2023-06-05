/* -*- c++ -*- ----------------------------------------------------------
Bespoke fix to compute the updated with valeus via a fixed point iteration,
and associated omega_SPH values for use in SPH force expressions.

Thomas Campbell
------------------------------------------------------------------------- */

#include "fix_dynamic_widths.h"
#include <cstring>
#include <mpi.h>
#include "group.h"
#include "modify.h"
#include "error.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "pair.h"
#include "pair_hybrid.h"
#include <cmath>
#include <math.h>
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "irregular.h"
#include "compute.h"
#include "domain.h"
#include "utils.h"
#include "memory.h"

using namespace LAMMPS_NS;
using namespace FixConst;

FixDynamicWidths::FixDynamicWidths(LAMMPS *lmp, int narg, char **arg):
	Fix(lmp, narg, arg)
{
  if (narg < 6) error->all(FLERR,"Illegal fix Dynamic Widths command");

  constant	= force->numeric(FLERR,arg[3]);
  N_iter = force->numeric(FLERR,arg[4]);
  mix_fact = force->numeric(FLERR,arg[5]);
  start_width = force->numeric(FLERR,arg[6]);
  cut_global = force->numeric(FLERR,arg[7]);
  pair_name = strdup(arg[8]);

  cutsquared = cut_global*cut_global;

  comm_forward = 1; 

}

FixDynamicWidths::~FixDynamicWidths()
{
  //free pair_name variable
  free(pair_name);
}

int FixDynamicWidths::setmask()
{
  int mask = 0;
  mask |= FixConst::POST_INTEGRATE;
  mask |= FixConst::POST_NEIGHBOR;
  mask |= FixConst::PRE_FORCE;
  return mask;
}

void FixDynamicWidths::init()
{ 
  double *width_SPH = atom->width_SPH;
  double *rho_SPH = atom->rho_SPH;
  double *omega_SPH = atom->omega_SPH;

  int natoms = static_cast<int> (atom->natoms);

  fprintf(screen,"FixDynamicWidths::init() being called... \n");
  fprintf(screen,"start_width =  %f\n", start_width);

  // assign all particles same initial width
  for(int i = 0; i <+ natoms; ++i){
    width_SPH[i] = start_width;
  }
}

void FixDynamicWidths::setup_post_neighbor()
{
  fprintf(screen,"FixDynamicWidths::setup_post_neighbor() being called... \n");

  // inherit neighbour lists from pair style
  pair = force->pair;
  // If a hybrid style is used we need to acces the correct sub-style.
  PairHybrid *hybrid_pair = dynamic_cast<PairHybrid*> (pair);

  if (hybrid_pair) {
    // The pair style is a hybrid style.
    if (!pair_name) error->all(FLERR,"When a hybrid pair-style is used, 'pair_name' must be set for the lagrangian solver.");
    int nstyles = hybrid_pair->nstyles;
    int found = -1;
    for (int m = 0; m < nstyles; ++m) {
      if (strcmp(pair_name, hybrid_pair->keywords[m]) == 0) {
        if (found != -1) error->all(FLERR,"Multiple pair-styles with the name 'pair_name' found.");
        found = m;
      }
    }
    if (found == -1) error->all(FLERR,"No pair-style with the name of 'pair_name' found.");
    // Set the correct pair style.
    pair = hybrid_pair->styles[found];
  }
}

void FixDynamicWidths::setup_pre_force(int)
{
  int a,i,j,ii,jj,inum,jnum,itype,jtype;
  int *ilist,*jlist,*numneigh,**firstneigh;
	double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *rho_SPH = atom->rho_SPH;
  double *omega_SPH = atom->omega_SPH;
  double *width_SPH = atom->width_SPH;
  double *mass = atom->mass;

  double pi_fact;
  double imass,jmass;
  double xtmp,ytmp,ztmp;
  double h_i,hm2_i,h_j,hm2_j;
  double delx,dely,delz;
  double gauss_pre_i,gauss_pre_j;
  double rsq;
  double dh_drho_SPH_i,dh_drho_SPH_j;
  double m_gauss_ij,m_gauss_ji;


	int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int newton_pair = force->newton_pair;

  fprintf(screen,"FixDynamicWidths::setup_pre_force() being called... \n");

  list = pair->list;
	
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  pi_fact = 1/pow(2*M_PI,1.5);

  int natoms = static_cast<int> (atom->natoms);

  // loop over N iterations as requested:

  for (a = 0; a < N_iter; a++){

    // clear all density values

    for(int i = 0; i <+ natoms; ++i){
      rho_SPH[i] = 0.;
    }

    // loop over my atoms

    for (ii = 0; ii < inum; ii++) {

      i = ilist[ii];

      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];

      itype = type[i];
      jlist = firstneigh[i];

      jnum = numneigh[i];
      imass = mass[itype];
      h_i = width_SPH[i];
      hm2_i = 1/(h_i*h_i);

      // 3D Gaussian prefactor
      gauss_pre_i = pi_fact*(1/(h_i*h_i*h_i));

      // self interaction term
      rho_SPH[i] += imass*gauss_pre_i;

      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;

        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];

        rsq = delx*delx + dely*dely + delz*delz;

        jtype = type[j];
  
        if (rsq < cutsquared) {

          jmass = mass[jtype];
          m_gauss_ij = jmass*gauss_pre_i*exp(-(rsq)*hm2_i/2);
          rho_SPH[i] += m_gauss_ij;

          if (newton_pair || j < nlocal) {


            h_j = width_SPH[j];
            hm2_j = 1/(h_j*h_j);
            gauss_pre_j = pi_fact*(1/(h_j*h_j*h_j));
            m_gauss_ji = imass*gauss_pre_j*exp(-(rsq)*hm2_j/2);
            rho_SPH[j] += m_gauss_ji;

          }
        }
      }
    }

    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      itype = type[i];
      imass = mass[itype];
      // mixing factor applied

      fprintf(screen,"imass =  %f\n", imass);
      fprintf(screen,"rho_SPH[i] =  %f\n", rho_SPH[i]);
      fprintf(screen,"pow(imass/rho_SPH[i],(1./3.)) =  %f\n", pow(imass/rho_SPH[i],(1./3.)));
      
      width_SPH[i] = mix_fact*constant*(pow(imass/rho_SPH[i],(1./3.))) + (1-mix_fact)*(width_SPH[i]);
      fprintf(screen,"width_SPH[i] =  %f\n", width_SPH[i]);
    }
  }

  // clear density (again) and omega_SPH values

  for(int i = 0; i <+ natoms; ++i){
      rho_SPH[i] = 0.;
      omega_SPH[i] = 0.;
    }

  // final loop to compute density values with mix factor widths
  
  for (ii = 0; ii < inum; ii++) {

    i = ilist[ii];

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    imass = mass[itype];

    h_i = width_SPH[i];
    hm2_i = 1/(h_i*h_i);

    // 3D Gaussian prefactor
    gauss_pre_i = pi_fact*(1/(h_i*h_i*h_i));

    // self interaction term
    rho_SPH[i] += imass*gauss_pre_i;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];

      rsq = delx*delx + dely*dely + delz*delz;

      jtype = type[j];

      if (rsq < cutsquared) {

        jmass = mass[jtype];
        m_gauss_ij = jmass*gauss_pre_i*exp(-(rsq)*hm2_i/2);
        rho_SPH[i] += m_gauss_ij;

        if (newton_pair || j < nlocal) {

          h_j = width_SPH[j];
          hm2_j = 1/(h_j*h_j);
          gauss_pre_j = pi_fact*(1/(h_j*h_j*h_j));
          m_gauss_ji = imass*gauss_pre_j*exp(-(rsq)*hm2_j/2);
          rho_SPH[j] += m_gauss_ji;

        }
      }
    }
  }
  // loop to compute omega_SPH values after assignment of density and width values
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    itype = type[i];
    jlist = firstneigh[i];

    jnum = numneigh[i];

    imass = mass[itype];

    h_i = width_SPH[i];

    // 3D Gaussian prefactor
    gauss_pre_i = pi_fact*(1/(h_i*h_i*h_i));
    // fprintf(screen,"i =  %d\n", i);
    // fprintf(screen,"rho_SPH[i] =  %f\n", rho_SPH[i]);


    dh_drho_SPH_i = -width_SPH[i]/(3*rho_SPH[i]);

    // 1 and self interaction term
    omega_SPH[i] = 1 - dh_drho_SPH_i*imass*(-3*gauss_pre_i/h_i);

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];

      rsq = delx*delx + dely*dely + delz*delz;

      jtype = type[j];

      if (rsq < cutsquared) {
        fprintf(screen,"gauss_pre_i =  %.11f\n", gauss_pre_i);
        fprintf(screen,"h_i =  %.11f\n", h_i);
        fprintf(screen,"rsq =  %.11f\n", rsq);
        fprintf(screen,"Gauss_Width_Deriv(gauss_pre_i,h_i,rsq) =  %.11f\n", Gauss_Width_Deriv(gauss_pre_i,h_i,rsq));


        jmass = mass[jtype];
        omega_SPH[i] -= dh_drho_SPH_i*jmass*Gauss_Width_Deriv(gauss_pre_i,h_i,rsq);

        if (newton_pair || j < nlocal) {

          h_j = width_SPH[j];
          gauss_pre_j = pi_fact*(1/(h_j*h_j*h_j));
          dh_drho_SPH_j = -width_SPH[j]/(3*rho_SPH[j]);
          omega_SPH[j] -= dh_drho_SPH_j*imass*Gauss_Width_Deriv(gauss_pre_j,h_j,rsq);

        }
      }
    }
  }
}

void FixDynamicWidths::pre_force(int)
{
  return;
}

void FixDynamicWidths::post_neighbor()
{
  return;
}

void FixDynamicWidths::post_integrate()
{
  int a,i,j,ii,jj,inum,jnum,itype,jtype;
  int *ilist,*jlist,*numneigh,**firstneigh;
	double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *rho_SPH = atom->rho_SPH;
  double *omega_SPH = atom->omega_SPH;
  double *width_SPH = atom->width_SPH;
  double *mass = atom->mass;

  double pi_fact;
  double imass,jmass;
  double xtmp,ytmp,ztmp;
  double h_i,hm2_i,h_j,hm2_j;
  double delx,dely,delz;
  double gauss_pre_i,gauss_pre_j;
  double rsq;
  double dh_drho_SPH_i,dh_drho_SPH_j;
  double m_gauss_ij,m_gauss_ji;


	int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int newton_pair = force->newton_pair;

  fprintf(screen,"FixDynamicWidths::post_integrate() being called check... \n");

  list = pair->list;
	
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  pi_fact = 1/pow(2*M_PI,1.5);

  int natoms = static_cast<int> (atom->natoms);

  // loop over N iterations as requested:

  for (a = 0; a < N_iter; a++){

    // clear all density values

    for(int i = 0; i <+ natoms; ++i){
      rho_SPH[i] = 0.;
    }

    // loop over my atoms

    for (ii = 0; ii < inum; ii++) {

      i = ilist[ii];

      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];

      itype = type[i];
      jlist = firstneigh[i];

      jnum = numneigh[i];
      imass = mass[itype];
      h_i = width_SPH[i];
      hm2_i = 1/(h_i*h_i);

      // 3D Gaussian prefactor
      gauss_pre_i = pi_fact*(1/(h_i*h_i*h_i));

      // self interaction term
      rho_SPH[i] += imass*gauss_pre_i;

      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;

        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];

        rsq = delx*delx + dely*dely + delz*delz;

        jtype = type[j];
  
        if (rsq < cutsquared) {

          jmass = mass[jtype];
          m_gauss_ij = jmass*gauss_pre_i*exp(-(rsq)*hm2_i/2);
          rho_SPH[i] += m_gauss_ij;

          if (newton_pair || j < nlocal) {


            h_j = width_SPH[j];
            hm2_j = 1/(h_j*h_j);
            gauss_pre_j = pi_fact*(1/(h_j*h_j*h_j));
            m_gauss_ji = imass*gauss_pre_j*exp(-(rsq)*hm2_j/2);
            rho_SPH[j] += m_gauss_ji;

          }
        }
      }
    }

    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      itype = type[i];
      imass = mass[itype];
      // mixing factor applied

      // fprintf(screen,"imass =  %f\n", imass);
      // fprintf(screen,"rho_SPH[i] =  %f\n", rho_SPH[i]);
      // fprintf(screen,"pow(imass/rho_SPH[i],(1./3.)) =  %f\n", pow(imass/rho_SPH[i],(1./3.)));
      
      width_SPH[i] = mix_fact*constant*(pow(imass/rho_SPH[i],(1./3.))) + (1-mix_fact)*(width_SPH[i]);
      // fprintf(screen,"width_SPH[i] =  %f\n", width_SPH[i]);
    }
  }

  // clear density (again) and omega_SPH values

  for(int i = 0; i <+ natoms; ++i){
      rho_SPH[i] = 0.;
      omega_SPH[i] = 0.;
    }

  // final loop to compute density values with mix factor widths
  
  for (ii = 0; ii < inum; ii++) {

    i = ilist[ii];

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    imass = mass[itype];

    h_i = width_SPH[i];
    hm2_i = 1/(h_i*h_i);

    // 3D Gaussian prefactor
    gauss_pre_i = pi_fact*(1/(h_i*h_i*h_i));

    // self interaction term
    rho_SPH[i] += imass*gauss_pre_i;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];

      rsq = delx*delx + dely*dely + delz*delz;

      jtype = type[j];

      if (rsq < cutsquared) {

        jmass = mass[jtype];
        m_gauss_ij = jmass*gauss_pre_i*exp(-(rsq)*hm2_i/2);
        rho_SPH[i] += m_gauss_ij;

        if (newton_pair || j < nlocal) {

          h_j = width_SPH[j];
          hm2_j = 1/(h_j*h_j);
          gauss_pre_j = pi_fact*(1/(h_j*h_j*h_j));
          m_gauss_ji = imass*gauss_pre_j*exp(-(rsq)*hm2_j/2);
          rho_SPH[j] += m_gauss_ji;

        }
      }
    }
  }
  comm->forward_comm_pair(this);
  // loop to compute omega_SPH values after assignment of density and width values
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    itype = type[i];
    jlist = firstneigh[i];

    jnum = numneigh[i];

    imass = mass[itype];

    h_i = width_SPH[i];

    // 3D Gaussian prefactor
    gauss_pre_i = pi_fact*(1/(h_i*h_i*h_i));
    // fprintf(screen,"i =  %d\n", i);
    // fprintf(screen,"rho_SPH[i] =  %f\n", rho_SPH[i]);

    dh_drho_SPH_i = -width_SPH[i]/(3*rho_SPH[i]);

    // 1 and self interaction term
    omega_SPH[i] = 1 - dh_drho_SPH_i*imass*(-3*gauss_pre_i/h_i);

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];

      rsq = delx*delx + dely*dely + delz*delz;


      jtype = type[j];

      if (rsq < cutsquared) {

        jmass = mass[jtype];

        // fprintf(screen,"gauss_pre_i =  %.11f\n", gauss_pre_i);
        // fprintf(screen,"h_i =  %.11f\n", h_i);
        // fprintf(screen,"rsq =  %.11f\n", rsq);
        // fprintf(screen,"Gauss_Width_Deriv(gauss_pre_i,h_i,rsq) =  %.11f\n", Gauss_Width_Deriv(gauss_pre_i,h_i,rsq));

        omega_SPH[i] -= dh_drho_SPH_i*jmass*Gauss_Width_Deriv(gauss_pre_i,h_i,rsq);

        if (newton_pair || j < nlocal) {

          h_j = width_SPH[j];
          gauss_pre_j = pi_fact*(1/(h_j*h_j*h_j));
          dh_drho_SPH_j = -width_SPH[j]/(3*rho_SPH[j]);
          omega_SPH[j] -= dh_drho_SPH_j*imass*Gauss_Width_Deriv(gauss_pre_j,h_j,rsq);

        }
      }
    }
  }
}
double FixDynamicWidths::Gauss_Width_Deriv(double pre_fact, double wid, double sep_sq)
{
  return pre_fact*exp(-sep_sq/(wid*wid*2.0))*((sep_sq)/(wid*wid*wid) - 3.0/wid);
}

int FixDynamicWidths::pack_forward_comm(int n, int *list, double *buf,
                               int /*pbc_flag*/, int * /*pbc*/)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = rho_SPH[j];
  }
    return m;
}

/* ---------------------------------------------------------------------- */

void FixDynamicWidths::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++){
    rho_SPH[i] = buf[m++];
  }
  
}

/* ----------------------------------------------------------------------
   init specific to this fix
------------------------------------------------------------------------- */

// void PairCoulCut::init_style()
// {
//   if (!atom->TC_SPH_flag)
//     error->all(FLERR,"fix_dynamic_widths requires atom attributes rho_SPH, width_SPH");

//   neighbor->request(this,instance_me);
// }