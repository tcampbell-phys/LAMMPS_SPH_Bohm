/* -*- c++ -*- ----------------------------------------------------------
Bespoke fix to compute the updated with valeus via a fixed point iteration,
and associated omega_SPH values for use in SPH force expressions.

Thomas Campbell
------------------------------------------------------------------------- */

#include "fix_fixed_widths.h"
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

FixFixedWidths::FixFixedWidths(LAMMPS *lmp, int narg, char **arg):
	Fix(lmp, narg, arg)
{
  if (narg < 6) error->all(FLERR,"Illegal fix Fixed Widths command");

  global_width = force->numeric(FLERR,arg[3]);
  cut_global = force->numeric(FLERR,arg[4]);
  // type not involved in SPH loop
  type_avoid = force->numeric(FLERR,arg[5]);
  pair_name = strdup(arg[6]);

  cutsquared = cut_global*cut_global;

  comm_forward = 1;
  comm_reverse = 1; 

  commflag = 0;

  assignment_flag = 1;

}

FixFixedWidths::~FixFixedWidths()
{
  //free pair_name variable
  free(pair_name);
}

int FixFixedWidths::setmask()
{
  int mask = 0;
  mask |= FixConst::POST_NEIGHBOR;
  mask |= FixConst::PRE_FORCE;
  mask |= FixConst::MIN_PRE_FORCE;
  return mask;
}

void FixFixedWidths::init()
{ 
  double *width_SPH = atom->width_SPH;
  double *rho_SPH = atom->rho_SPH;
  double *omega_SPH = atom->omega_SPH;
  double *u_SPH = atom->u_SPH;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  // assign all particles same initial width if initialising for the first time
  if (assignment_flag == 1){
    for(int i = 0; i < nall; ++i){
      width_SPH[i] = global_width;
      omega_SPH[i] = 1.;
      u_SPH[i] = 0.0;
    }
  }
  assignment_flag = 0;

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


void FixFixedWidths::setup_pre_force(int)
{
  FixedPointIterator();
}

void FixFixedWidths::min_pre_force(int)
{
  FixedPointIterator();
}


void FixFixedWidths::setup_post_neighbor()
{
  // inherit neighbour lists from pair style
  
  FixedPointIterator();
}

void FixFixedWidths::post_neighbor()
{
  return;
}

void FixFixedWidths::pre_force(int)
{
  FixedPointIterator();
}

void FixFixedWidths::FixedPointIterator()
{
  int a,i,j,ii,jj,inum,jnum,itype,jtype;
  int *ilist,*jlist,*numneigh,**firstneigh;
	double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *rho_SPH = atom->rho_SPH;
  double *dx_rho_SPH = atom->dx_rho_SPH;
  double *dy_rho_SPH = atom->dy_rho_SPH;
  double *dz_rho_SPH = atom->dz_rho_SPH;
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

  list = pair->list;
	
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  pi_fact = 1/pow(2*M_PI,1.5);

  // clear all density values

  for(int i = 0; i < nall; i++){
    rho_SPH[i] = 0.;
    dx_rho_SPH[i] = 0.;
    dy_rho_SPH[i] = 0.;
    dz_rho_SPH[i] = 0.;
    omega_SPH[i] = 1.;
    width_SPH[i] = global_width;
  }

  // loop to compute density values

  h_i = global_width;
  hm2_i = 1/(h_i*h_i);

  // 3D Gaussian prefactor
  gauss_pre_i = pi_fact*(1/(h_i*h_i*h_i));

  h_j = global_width;
  hm2_j = 1/(h_j*h_j);

  gauss_pre_j = pi_fact*(1/(h_j*h_j*h_j));
  
  for (ii = 0; ii < inum; ii++) {

    i = ilist[ii];

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    itype = type[i];
    if (type[i] == type_avoid){
      continue;
    }
    
    jlist = firstneigh[i];
    jnum = numneigh[i];
    imass = mass[itype];

    // self interaction term
    rho_SPH[i] += imass*gauss_pre_i;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      jtype = type[j];
      if (type[j] == type_avoid){
        continue;
      }

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];

      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < cutsquared) {

        jmass = mass[jtype];
        m_gauss_ij = jmass*gauss_pre_i*exp(-(rsq)*hm2_i/2);
        rho_SPH[i] += m_gauss_ij;
        dx_rho_SPH[i] += ((-delx)*hm2_i)*m_gauss_ij;
        dy_rho_SPH[i] += ((-dely)*hm2_i)*m_gauss_ij;
        dz_rho_SPH[i] += ((-delz)*hm2_i)*m_gauss_ij;

        if (newton_pair || j < nlocal) {

          m_gauss_ji = imass*gauss_pre_j*exp(-(rsq)*hm2_j/2);
          rho_SPH[j] += m_gauss_ji;
          dx_rho_SPH[j] += ((delx)*hm2_j)*m_gauss_ji;
          dy_rho_SPH[j] += ((dely)*hm2_j)*m_gauss_ji;
          dz_rho_SPH[j] += ((delz)*hm2_j)*m_gauss_ji;
          
        }
      }
    }
  }
  commflag = 3;
  comm_forward = 4;
  comm_reverse = 4; 
  if (newton_pair){
    comm->reverse_comm_fix(this);
  }
  comm->forward_comm_fix(this);
}

double FixFixedWidths::Gauss_Width_Deriv(double pre_fact, double wid, double sep_sq)
{
  return pre_fact*exp(-sep_sq/(wid*wid*2.0))*((sep_sq)/(wid*wid*wid) - 3.0/wid);
}

int FixFixedWidths::pack_forward_comm(int n, int *list, double *buf,
                               int /*pbc_flag*/, int * /*pbc*/)
{
  int i,j,m;

  m = 0;
  if (commflag == 3){
    double *rho_SPH = atom->rho_SPH;
    double *dx_rho_SPH = atom->dx_rho_SPH;
    double *dy_rho_SPH = atom->dy_rho_SPH;
    double *dz_rho_SPH = atom->dz_rho_SPH;
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = rho_SPH[j];
      buf[m++] = dx_rho_SPH[j];
      buf[m++] = dy_rho_SPH[j];
      buf[m++] = dz_rho_SPH[j];
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixFixedWidths::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  if (commflag == 3){
    double *rho_SPH = atom->rho_SPH;
    double *dx_rho_SPH = atom->dx_rho_SPH;
    double *dy_rho_SPH = atom->dy_rho_SPH;
    double *dz_rho_SPH = atom->dz_rho_SPH;
    for (i = first; i < last; i++){
      rho_SPH[i] = buf[m++];
      dx_rho_SPH[i] = buf[m++];
      dy_rho_SPH[i] = buf[m++];
      dz_rho_SPH[i] = buf[m++];
    }
  }
  
}

int FixFixedWidths::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  if (commflag == 3){
    double *rho_SPH = atom->rho_SPH;
    double *dx_rho_SPH = atom->dx_rho_SPH;
    double *dy_rho_SPH = atom->dy_rho_SPH;
    double *dz_rho_SPH = atom->dz_rho_SPH;
    for (i = first; i < last; i++){
      buf[m++] = rho_SPH[i];
      buf[m++] = dx_rho_SPH[i];
      buf[m++] = dy_rho_SPH[i];
      buf[m++] = dz_rho_SPH[i];
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixFixedWidths::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  if (commflag == 3){
    double *rho_SPH = atom->rho_SPH;
    double *dx_rho_SPH = atom->dx_rho_SPH;
    double *dy_rho_SPH = atom->dy_rho_SPH;
    double *dz_rho_SPH = atom->dz_rho_SPH;
    for (i = 0; i < n; i++) {
      j = list[i];
      rho_SPH[j] += buf[m++];
      dx_rho_SPH[j] += buf[m++];
      dy_rho_SPH[j] += buf[m++];
      dz_rho_SPH[j] += buf[m++];
    }
  }
}