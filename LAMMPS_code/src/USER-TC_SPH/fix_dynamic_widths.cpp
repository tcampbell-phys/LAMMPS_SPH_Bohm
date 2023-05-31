/* -*- c++ -*- ----------------------------------------------------------
Bespoke fix to compute the updated with valeus via a fixed point iteration,
and associated omega_SPH values for use in SPH force expressions.

Thomas Campbell
------------------------------------------------------------------------- */

#include "fix_dynamic_widths.h"
#include <cstring>
#include "group.h"
#include "modify.h"
#include "error.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include <cmath>
#include <math.h>
#include "comm.h"
#include "neighbor.h"
#include "irregular.h"
#include "fix_deform.h"
#include "compute.h"
#include "domain.h"
#include "memory.h"

using namespace LAMMPS_NS;
using namespace FixConst;

FixDynamicWidths::FixDynamicWidths(LAMMPS *lmp, int narg, char **arg):
	Fix(lmp, narg, arg), id_temp(NULL)
{
  if (narg < 5) error->all(FLERR,"Illegal fix Dynamic Widths command");

  constant	= force->numeric(FLERR,arg[3]);
  N_iter = force->numeric(FLERR,arg[4]);
  mix_fact = force->numeric(FLERR,arg[5]);
  start_width = force->numeric(FLERR,arg[6]);
  cut_global = force->numeric(FLERR,arg[7]);

  cutsquared = cut_global*cut_global; 

}	

int FixDynamicWidths::setmask()
{
  int mask = 0;
  mask |= FixConst::POST_INTEGRATE;
  return mask;
}

void FixDynamicWidths::init()
{
  double *width_SPH = atom->width_SPH;

  int natoms = static_cast<int> (atom->natoms);
  
  // assign all particles same initial width
  for(int i = 0; i <+ natoms; ++i){
    width_SPH[i] = start_width;
  }
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
  double dh_drho_SPH_i;
  double m_gauss_ij,m_gauss_ji;


	int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int newton_pair = force->newton_pair;
	
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
          h_j = width_SPH[j];
          hm2_j = 1/(h_j*h_j);
          gauss_pre_j = pi_fact*(1/(h_j*h_j*h_j));
          m_gauss_ij = jmass*gauss_pre_i*exp(-(rsq)*hm2_i/2);
          rho_SPH[i] += m_gauss_ij;

          if (newton_pair || j < nlocal) {

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
      width_SPH[i] = mix_fact*constant*(pow(imass/rho_SPH[i],(1/3))) + (1-mix_fact)*(width_SPH[i]);
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
        h_j = width_SPH[j];
        hm2_j = 1/(h_j*h_j);
        gauss_pre_j = pi_fact*(1/(h_j*h_j*h_j));
        m_gauss_ij = jmass*gauss_pre_i*exp(-(rsq)*hm2_i/2);
        rho_SPH[i] += m_gauss_ij;

        if (newton_pair || j < nlocal) {

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

    dh_drho_SPH_i = -width_SPH[i]/(3*rho_SPH[i]);

    // 1 and self interaction term
    omega_SPH[i] = 1 - dh_drho_SPH_i*(-3*gauss_pre_i/h_i);

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
        omega_SPH[i] -= dh_drho_SPH_i*jmass*Gauss_Width_Deriv(gauss_pre_i,h_i,rsq)

      }
    }
  }
}

double FixDynamicWidths::Gauss_Width_Deriv(double pre_fact, double wid, double sep_sq)
{
  return pre_fact*exp(-sep_sq/(wid*wid*2))*((sep_sq)/(wid*wid*wid) - 3/wid);
}