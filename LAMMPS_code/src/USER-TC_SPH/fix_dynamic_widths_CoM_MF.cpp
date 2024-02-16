/* -*- c++ -*- ----------------------------------------------------------
Bespoke fix to compute the updated with valeus via a fixed point iteration,
and associated omega_SPH values for use in SPH force expressions.

Widths updated according to local density of member electron. (Compatible with MF Bohm).

Thomas Campbell
------------------------------------------------------------------------- */

#include "fix_dynamic_widths_CoM_MF.h"
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

FixDynamicWidthsCoMMF::FixDynamicWidthsCoMMF(LAMMPS *lmp, int narg, char **arg):
	Fix(lmp, narg, arg), mu(NULL), mu_all(NULL), tau(NULL), tau_all(NULL)
{
  if (narg < 13) error->all(FLERR,"Illegal fix Dynamic Widths command");

  constant	= force->numeric(FLERR,arg[3]);
  N_iter = force->numeric(FLERR,arg[4]);
  mix_fact = force->numeric(FLERR,arg[5]);
  start_width = force->numeric(FLERR,arg[6]);
  cut_global = force->numeric(FLERR,arg[7]);
  // type not involved in SPH loop
  box_len = force->numeric(FLERR,arg[8]); //Accurate for a cubic structure
  type_avoid = force->numeric(FLERR,arg[9]);
  N_electrons = force->numeric(FLERR,arg[10]);
  N_elements_per_electron = force->numeric(FLERR,arg[11]);
  tag_ele_start = force->numeric(FLERR,arg[12]);

  pair_name = strdup(arg[13]);

  cutsquared = cut_global*cut_global;

  two_pi_over_len = 2*M_PI/box_len;

  // x_mu = NULL;
  // x_tau = NULL;
  // y_mu = NULL;
  // y_tau = NULL;
  // z_mu = NULL;
  // z_tau = NULL;

  mu_all = NULL;
  tau_all = NULL;

  unallocated = 1;
  allocated = 0;

  comm_forward = 1;
  comm_reverse = 1; 

  commflag = 0;
  nmax = 0;

  assignment_flag = 1;
  fprintf(screen,"\nIn FixDynamicWidthsCoMMF with MPI_Allreduce and no commflag = 4 or 5");
}

FixDynamicWidthsCoMMF::~FixDynamicWidthsCoMMF()
{
  //free pair_name variable
  free(pair_name);
  if (allocated){
    deallocate();
  }
}

int FixDynamicWidthsCoMMF::setmask()
{
  int mask = 0;
  mask |= FixConst::POST_NEIGHBOR;
  mask |= FixConst::PRE_FORCE;
  mask |= FixConst::MIN_PRE_FORCE;
  return mask;
}

void FixDynamicWidthsCoMMF::init()
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
      width_SPH[i] = start_width;
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


void FixDynamicWidthsCoMMF::setup_pre_force(int)
{
  FixedPointIterator();
  CoM_Calculator();
}

void FixDynamicWidthsCoMMF::min_pre_force(int)
{
  FixedPointIterator();
  CoM_Calculator();
}


void FixDynamicWidthsCoMMF::setup_post_neighbor()
{
  // inherit neighbour lists from pair style
  FixedPointIterator();
  CoM_Calculator();
}

void FixDynamicWidthsCoMMF::post_neighbor()
{
  return;
}

void FixDynamicWidthsCoMMF::pre_force(int)
{
  FixedPointIterator();
  CoM_Calculator();
}

void FixDynamicWidthsCoMMF::FixedPointIterator()
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

  int *tagid = atom->tag;
  int lo_lim_lev;
  int hi_lim_lev;

  list = pair->list;
	
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  pi_fact = 1/pow(2*M_PI,1.5);

  // loop over N iterations as requested:

  for (a = 0; a < N_iter; a++){

    // clear all density values

    for(int b = 0; b < nall; b++){
      rho_SPH[b] = 0.;
    }

    // loop over my atoms

    for (ii = 0; ii < inum; ii++) {

      i = ilist[ii];

      itype = type[i];
      if (type[i] == type_avoid){
        continue;
      }

      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];

      // determine owner electron ID limits
      lo_lim_lev = floor((tagid[i]-tag_ele_start)/N_elements_per_electron)*N_elements_per_electron + tag_ele_start;
      hi_lim_lev = lo_lim_lev + N_elements_per_electron;
      
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

        jtype = type[j];
        if (type[j] == type_avoid){
          continue;
        }

        // skip elements not in owner electron
        if ( tagid[j] >= hi_lim_lev || tagid[j] < lo_lim_lev ){
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
    commflag = 0;
    if (newton_pair){
      comm->reverse_comm_fix(this);
    }
    comm->forward_comm_fix(this);

    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      itype = type[i];
      if (type[i] == type_avoid){
        continue;
      }
      
      imass = mass[itype];
      // mixing factor applied
      width_SPH[i] = mix_fact*constant*(pow(imass/rho_SPH[i],(1./3.))) + (1-mix_fact)*(width_SPH[i]);
    }
    commflag = 1;
    comm->forward_comm_fix(this);
  }

  // clear density (again)

  for(int i = 0; i < nall; i++){
      
      rho_SPH[i] = 0.;
      dx_rho_SPH[i] = 0.;
      dy_rho_SPH[i] = 0.;
      dz_rho_SPH[i] = 0.;
      omega_SPH[i] = 0.;
    }

  // final loop to compute density values with mix factor widths
  
  for (ii = 0; ii < inum; ii++) {

    i = ilist[ii];

    itype = type[i];
    if (type[i] == type_avoid){
      continue;
    }

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    // determine owner electron ID limits
    lo_lim_lev = floor((tagid[i]-tag_ele_start)/N_elements_per_electron)*N_elements_per_electron + tag_ele_start;
    hi_lim_lev = lo_lim_lev + N_elements_per_electron;
    
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

      jtype = type[j];
      if (type[j] == type_avoid){
        continue;
      }

      // skip elements not in owner electron
      if ( tagid[j] >= hi_lim_lev || tagid[j] < lo_lim_lev ){
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

          h_j = width_SPH[j];
          hm2_j = 1/(h_j*h_j);
          gauss_pre_j = pi_fact*(1/(h_j*h_j*h_j));
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

  // loop to compute omega_SPH values after assignment of density and width values
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];

    // // // fprintf(screen,"\nrho_SPH[%d] = %16.16f",i,rho_SPH[i]);

    itype = type[i];
    if (type[i] == type_avoid){
      continue;
    }
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    // determine owner electron ID limits
    lo_lim_lev = floor((tagid[i]-tag_ele_start)/N_elements_per_electron)*N_elements_per_electron + tag_ele_start;
    hi_lim_lev = lo_lim_lev + N_elements_per_electron;
    
    jlist = firstneigh[i];

    jnum = numneigh[i];

    imass = mass[itype];

    h_i = width_SPH[i];

    // 3D Gaussian prefactor
    gauss_pre_i = pi_fact*(1/(h_i*h_i*h_i));

    dh_drho_SPH_i = -width_SPH[i]/(3*(rho_SPH[i]));

    // self interaction term
    omega_SPH[i] += 1. - dh_drho_SPH_i*imass*(-3*gauss_pre_i/h_i);

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      jtype = type[j];
      if (type[j] == type_avoid){
        continue;
      }

      // skip elements not in owner electron
      if ( tagid[j] >= hi_lim_lev || tagid[j] < lo_lim_lev ){
        continue;
      }

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];

      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < cutsquared) {

        jmass = mass[jtype];
        omega_SPH[i] -= dh_drho_SPH_i*jmass*Gauss_Width_Deriv(gauss_pre_i,h_i,rsq);

        if (newton_pair || j < nlocal) {

          h_j = width_SPH[j];
          gauss_pre_j = pi_fact*(1/(h_j*h_j*h_j));
          dh_drho_SPH_j = -width_SPH[j]/(3*(rho_SPH[j]));
          omega_SPH[j] -= dh_drho_SPH_j*imass*Gauss_Width_Deriv(gauss_pre_j,h_j,rsq);
        }
      }
    }
  }
  commflag = 2;
  comm_forward = 1;
  comm_reverse = 1; 
  if (newton_pair){
    comm->reverse_comm_fix(this);
  }
  comm->forward_comm_fix(this);
}

double FixDynamicWidthsCoMMF::Gauss_Width_Deriv(double pre_fact, double wid, double sep_sq)
{
  return pre_fact*exp(-sep_sq/(wid*wid*2.0))*((sep_sq)/(wid*wid*wid) - 3.0/wid);
}

void FixDynamicWidthsCoMMF::CoM_Calculator()
{
  int i,ii,inum;
  int *ilist;
	double **x = atom->x;
  double *x_COM = atom->x_COM;
  double *y_COM = atom->y_COM;
  double *z_COM = atom->z_COM;
  double *width_SPH = atom->width_SPH;
  double *mass = atom->mass;

  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  int *tagid = atom->tag;

  list = pair->list;

  // int ntimestep = update->ntimestep;
	
  inum = list->inum;
  ilist = list->ilist;

  double theta_ave;

  int ele_ind;

  if (unallocated){
    deallocate();
    // fprintf(screen,"\nFixDynamicWidthsCoMMF allocating mu tau memory...\n");
    allocate();
    allocated = 1;
    unallocated = 0;
  }

  for (int a = 0; a < 3*N_electrons; a++){ 
    mu[a] = 0.0;
    tau[a] = 0.0;
  }

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (type[i] == type_avoid){
      continue;
    }
    ele_ind = floor((tagid[i]-tag_ele_start)/N_elements_per_electron);
    mu[ele_ind] += cos(x[i][0]*two_pi_over_len);
    mu[ele_ind+N_electrons] += cos(x[i][1]*two_pi_over_len);
    mu[ele_ind+2*N_electrons] += cos(x[i][2]*two_pi_over_len);
    tau[ele_ind] += sin(x[i][0]*two_pi_over_len);
    tau[ele_ind+N_electrons] += sin(x[i][1]*two_pi_over_len);
    tau[ele_ind+2*N_electrons] += sin(x[i][2]*two_pi_over_len);
  }

  MPI_Allreduce(mu,mu_all,3*N_electrons,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(tau,tau_all,3*N_electrons,MPI_DOUBLE,MPI_SUM,world);

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (type[i] == type_avoid){
      continue;
    }
    ele_ind = floor((tagid[i]-tag_ele_start)/N_elements_per_electron);
    x_COM[i] = (atan2(-tau_all[ele_ind]/N_elements_per_electron,-mu_all[ele_ind]/N_elements_per_electron) + M_PI)/two_pi_over_len;
    y_COM[i] = (atan2(-tau_all[ele_ind+N_electrons]/N_elements_per_electron,-mu_all[ele_ind+N_electrons]/N_elements_per_electron) + M_PI)/two_pi_over_len;
    z_COM[i] = (atan2(-tau_all[ele_ind+2*N_electrons]/N_elements_per_electron,-mu_all[ele_ind+2*N_electrons]/N_elements_per_electron) + M_PI)/two_pi_over_len;
    // if (tagid[i] == 801){
    //   fprintf(screen,"\n\ntimestep %d \n\nx_COM = %16.16f \ny_COM = %16.16f \nz_COM = %16.16f",ntimestep,x_COM[i],y_COM[i],z_COM[i]);
    // }
  }
}

int FixDynamicWidthsCoMMF::pack_forward_comm(int n, int *list, double *buf,
                               int /*pbc_flag*/, int * /*pbc*/)
{
  int i,j,m;

  m = 0;
  if (commflag == 0){
    double *rho_SPH = atom->rho_SPH;
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = rho_SPH[j];
    }
  }
  if (commflag == 5){
    double *x_COM = atom->x_COM;
    double *y_COM = atom->y_COM;
    double *z_COM = atom->z_COM;
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x_COM[j];
      buf[m++] = y_COM[j];
      buf[m++] = z_COM[j];
      
    }
  }
  if (commflag == 4){
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x_mu[j];
      buf[m++] = x_tau[j];
      buf[m++] = y_mu[j];
      buf[m++] = y_tau[j];
      buf[m++] = z_mu[j];
      buf[m++] = z_tau[j];
    }
  }
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
  if (commflag == 1){
    double *width_SPH = atom->width_SPH;
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = width_SPH[j];
    }
  }
  if (commflag == 2){
    double *omega_SPH = atom->omega_SPH;
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = omega_SPH[j];
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixDynamicWidthsCoMMF::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;

  if (commflag == 5){
    double *x_COM = atom->x_COM;
    double *y_COM = atom->y_COM;
    double *z_COM = atom->z_COM;
    for (i = first; i < last; i++){
      x_COM[i] = buf[m++];
      y_COM[i] = buf[m++];
      z_COM[i] = buf[m++];
    }
  }

  if (commflag == 4){
    for (i = first; i < last; i++){
      x_mu[i] = buf[m++];
      x_tau[i] = buf[m++];
      y_mu[i] = buf[m++];
      y_tau[i] = buf[m++];
      z_mu[i] = buf[m++];
      z_tau[i] = buf[m++];
    }
  }
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
  if (commflag == 0){
    double *rho_SPH = atom->rho_SPH;
    for (i = first; i < last; i++){
      rho_SPH[i] = buf[m++];
    }
  }
  if (commflag == 1){
    double *width_SPH = atom->width_SPH;
    for (i = first; i < last; i++){
      width_SPH[i] = buf[m++];
    }
  }
  if (commflag == 2){
    double *omega_SPH = atom->omega_SPH;
    for (i = first; i < last; i++){
      omega_SPH[i] = buf[m++];
    }
  }
  
}

int FixDynamicWidthsCoMMF::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  if (commflag == 5){
    double *x_COM = atom->x_COM;
    double *y_COM = atom->y_COM;
    double *z_COM = atom->z_COM;
    for (i = first; i < last; i++){
      buf[m++] = x_COM[i];
      buf[m++] = y_COM[i];
      buf[m++] = z_COM[i];
    }
  }
  if (commflag == 4){
    for (i = first; i < last; i++){
      buf[m++] = x_mu[i];
      buf[m++] = x_tau[i];
      buf[m++] = y_mu[i];
      buf[m++] = y_tau[i];
      buf[m++] = z_mu[i];
      buf[m++] = z_tau[i];
    }
  }
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
  if (commflag == 0){
    double *rho_SPH = atom->rho_SPH;
    for (i = first; i < last; i++){
      buf[m++] = rho_SPH[i];
    }
  }
  if (commflag == 1){
    double *width_SPH = atom->width_SPH;
    for (i = first; i < last; i++){
      buf[m++] = width_SPH[i];
    }
  }
  if (commflag == 2){
    double *omega_SPH = atom->omega_SPH;
    for (i = first; i < last; i++){
      buf[m++] = omega_SPH[i];
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixDynamicWidthsCoMMF::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  if (commflag == 5){
    double *x_COM = atom->x_COM;
    double *y_COM = atom->y_COM;
    double *z_COM = atom->z_COM;
    for (i = 0; i < n; i++) {
      j = list[i];
      x_COM[j] += buf[m++];
      y_COM[j] += buf[m++];
      z_COM[j] += buf[m++];
    }
  }
  if (commflag == 4){
    for (i = 0; i < n; i++) {
      j = list[i];
      x_mu[j] += buf[m++];
      x_tau[j] += buf[m++];
      y_mu[j] += buf[m++];
      y_tau[j] += buf[m++];
      z_mu[j] += buf[m++];
      z_tau[j] += buf[m++];
    }
  }
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
  if (commflag == 0){
    double *rho_SPH = atom->rho_SPH;
    for (i = 0; i < n; i++) {
      j = list[i];
      rho_SPH[j] += buf[m++];
    }
  }
  if (commflag == 1){
    double *width_SPH = atom->width_SPH;
    for (i = 0; i < n; i++) {
      j = list[i];
      width_SPH[j] += buf[m++];
    }
  }
  if (commflag == 2){
    double *omega_SPH = atom->omega_SPH;
    for (i = 0; i < n; i++) {
      j = list[i];
      omega_SPH[j] += buf[m++];
    }
  }
}

/* ----------------------------------------------------------------------
   allocate memory that depends on # of electrons
------------------------------------------------------------------------- */

void FixDynamicWidthsCoMMF::allocate()
{
  fprintf(screen,"\nFixDynamicWidthsCoMMF::allocate()\n");
  // Centre of Mass array x 0:N_ele, y N_ele:2*N_ele, z 2*N_ele:3_Nele
  mu = new double[3*N_electrons]; 
  mu_all = new double[3*N_electrons];
  tau = new double[3*N_electrons]; 
  tau_all = new double[3*N_electrons];
}

/* ----------------------------------------------------------------------
   deallocate memory that depends on # of electrons
------------------------------------------------------------------------- */

void FixDynamicWidthsCoMMF::deallocate()
{
  fprintf(screen,"\nFixDynamicWidthsCoMMF::deallocate()\n");
  delete [] mu;
  delete [] mu_all;
  delete [] tau;
  delete [] tau_all;
}