
/* -*- c++ -*- ----------------------------------------------------------
Bespoke fix to compute the updated with valeus via a fixed point iteration,
and associated omega_SPH values for use in SPH force expressions.

Widths updated according to total density of electrons. (Compatible with QHD Bohm).

Thomas Campbell
------------------------------------------------------------------------- */

#include "fix_dynamic_widths_CoM.h"
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

FixDynamicWidthsCoMBPMoczRobFtwo::FixDynamicWidthsCoMBPMoczRobFtwo(LAMMPS *lmp, int narg, char **arg):
	Fix(lmp, narg, arg), mu(NULL), mu_all(NULL), tau(NULL), tau_all(NULL)
{
  if (narg < 14) error->all(FLERR,"Illegal fix Dynamic Widths command");

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
  gamma_factor = force->numeric(FLERR,arg[13]);

  pair_name = strdup(arg[14]);

  cutsquared = cut_global*cut_global;

  two_pi_over_len = 2*M_PI/box_len;

  nmax = 0;

  scalar_flag = 1;
  global_freq = 1;
  extscalar = 1;

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

  // manybody_flag = 1;

  dxx_rho_fix = NULL;
  dxy_rho_fix = NULL;
  dxz_rho_fix = NULL;
  dyy_rho_fix = NULL;
  dyz_rho_fix = NULL;
  dzz_rho_fix = NULL;

  // compute dx_rho_fix terms internally in Bohm algorithm as F2 derivatives

  dx_rho_fix = NULL;
  dy_rho_fix = NULL;
  dz_rho_fix = NULL;

  Pxx_fix = NULL;
  Pxy_fix = NULL;
  Pxz_fix = NULL;
  Pyz_fix = NULL;
  Pyy_fix = NULL;
  Pzz_fix = NULL;

  comm_forward = 9;
  comm_reverse = 9;

  hplanck  = force->hplanck;
  hbar = hplanck/(2*M_PI);
  // Bohm pressure prefactor
  f_prefactor = force->hhmrr2e * (hbar*hbar)/(4*e_mass);

  pi_fact = 1/pow(2*M_PI,1.5);
  // // fprintf(screen,"\nIn FixDynamicWidthsCoMBPMoczRobFtwo with MPI_Allreduce and no commflag = 4 or 5");

  eflag = 0;
  eBohm = 0.0;
}

FixDynamicWidthsCoMBPMoczRobFtwo::~FixDynamicWidthsCoMBPMoczRobFtwo()
{
  //free pair_name variable
  free(pair_name);
  if (allocated){
    deallocate();
    memory->destroy(dxx_rho_fix);
    memory->destroy(dxy_rho_fix);
    memory->destroy(dxz_rho_fix);
    memory->destroy(dyy_rho_fix);
    memory->destroy(dyz_rho_fix);
    memory->destroy(dzz_rho_fix);

    memory->destroy(dx_rho_fix);
    memory->destroy(dy_rho_fix);
    memory->destroy(dz_rho_fix);

    memory->destroy(Pxx_fix);
    memory->destroy(Pxy_fix);
    memory->destroy(Pxz_fix);
    memory->destroy(Pyy_fix);
    memory->destroy(Pyz_fix);
    memory->destroy(Pzz_fix);
  }
}

int FixDynamicWidthsCoMBPMoczRobFtwo::setmask()
{
  int mask = 0;
  mask |= FixConst::POST_NEIGHBOR;
  mask |= FixConst::PRE_FORCE;
  mask |= FixConst::MIN_PRE_FORCE;
  mask |= FixConst::END_OF_STEP;
  mask |= FixConst::THERMO_ENERGY;
  return mask;
}

void FixDynamicWidthsCoMBPMoczRobFtwo::init()
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


void FixDynamicWidthsCoMBPMoczRobFtwo::setup_pre_force(int)
{
  FixedPointIterator();
  CoM_Calculator();
}

void FixDynamicWidthsCoMBPMoczRobFtwo::min_pre_force(int)
{
  FixedPointIterator();
  CoM_Calculator();
}


void FixDynamicWidthsCoMBPMoczRobFtwo::setup_post_neighbor()
{
  // inherit neighbour lists from pair style
  FixedPointIterator();
  CoM_Calculator();
}

void FixDynamicWidthsCoMBPMoczRobFtwo::post_neighbor()
{
  return;
}

void FixDynamicWidthsCoMBPMoczRobFtwo::pre_force(int)
{
  FixedPointIterator();
  CoM_Calculator();
}

void FixDynamicWidthsCoMBPMoczRobFtwo::end_of_step()
{
  // fprintf(screen,"\nFixDynamicWidthsCoMBPMoczRobFtwo end_of_step");
  CalcBohmPot();
}

void FixDynamicWidthsCoMBPMoczRobFtwo::FixedPointIterator()
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

  // loop over N iterations as requested:

  for (a = 0; a < N_iter; a++){

    // clear all density values

    for(int b = 0; b < nall; b++){
      rho_SPH[b] = 0.;
    }

    // loop over my atoms

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

    // // // // // fprintf(screen,"\nrho_SPH[%d] = %16.16f",i,rho_SPH[i]);

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

double FixDynamicWidthsCoMBPMoczRobFtwo::Gauss_Width_Deriv(double pre_fact, double wid, double sep_sq)
{
  return pre_fact*exp(-sep_sq/(wid*wid*2.0))*((sep_sq)/(wid*wid*wid) - 3.0/wid);
}

void FixDynamicWidthsCoMBPMoczRobFtwo::CoM_Calculator()
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
    // // // fprintf(screen,"\nFixDynamicWidthsCoMBPMoczRobFtwo allocating mu tau memory...\n");
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
    //   // // fprintf(screen,"\n\ntimestep %d \n\nx_COM = %16.16f \ny_COM = %16.16f \nz_COM = %16.16f",ntimestep,x_COM[i],y_COM[i],z_COM[i]);
    // }
  }
}

void FixDynamicWidthsCoMBPMoczRobFtwo::CalcBohmPot()
{
  // fprintf(screen,"\nCheckpoint 0");
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz;
  double vxtmp,vytmp,vztmp,delvx,delvy,delvz;
  double xtmp2,ytmp2,ztmp2,delx2,dely2,delz2;
  double delx_2,dely_2,delz_2,rsq;
  double gauss_pre_i;
  double gauss_pre_j;
  double h_i,h2_i,hm2_i,hm4_i;
  double h_j,h2_j,hm2_j,hm4_j;
  double imass,jmass,ijmass;
  double m_gauss_ij,m_gauss_ji;
  double exp_ij,exp_ji;
  double omega_i,omega_j;
  double dx_Wij,dy_Wij,dz_Wij;
  double dx_Wji,dy_Wji,dz_Wji;
  double rho_i2,rho_j2;
  double dens_factor_ij,dens_factor_ji;
  double rho_i,rho_j;
  double bohm_pot;
  double u_prefact_i,u_prefact_j;

  // pair = force->pair;

  // double eng_vdwl = force->pair->eng_vdwl;
  // // fprintf(screen,"\neng_vdwl = %f",eng_vdwl);

  // ev_init(1,0);
  // fprintf(screen,"\nCheckpoint 1");

  eflag = 0;
  eBohm = 0.0;

  int *ilist,*jlist,*numneigh,**firstneigh;

  if (atom->nmax > nmax) {
    // delete and create new memory arrays for any per-particle variables that need communicating.
    memory->destroy(dxx_rho_fix);
    memory->destroy(dxy_rho_fix);
    memory->destroy(dxz_rho_fix);
    memory->destroy(dyy_rho_fix);
    memory->destroy(dyz_rho_fix);
    memory->destroy(dzz_rho_fix);

    memory->destroy(dx_rho_fix);
    memory->destroy(dy_rho_fix);
    memory->destroy(dz_rho_fix);

    memory->destroy(Pxx_fix);
    memory->destroy(Pxy_fix);
    memory->destroy(Pxz_fix);
    memory->destroy(Pyy_fix);
    memory->destroy(Pyz_fix);
    memory->destroy(Pzz_fix);
    
    nmax = atom->nmax;

    memory->create(dxx_rho_fix,nmax,"fix:dxx_rho_fix");
    memory->create(dxy_rho_fix,nmax,"fix:dxy_rho_fix");
    memory->create(dxz_rho_fix,nmax,"fix:dxz_rho_fix");
    memory->create(dyy_rho_fix,nmax,"fix:dyy_rho_fix");
    memory->create(dyz_rho_fix,nmax,"fix:dyz_rho_fix");
    memory->create(dzz_rho_fix,nmax,"fix:dzz_rho_fix");

    memory->create(dx_rho_fix,nmax,"fix:dx_rho_fix");
    memory->create(dy_rho_fix,nmax,"fix:dy_rho_fix");
    memory->create(dz_rho_fix,nmax,"fix:dz_rho_fix");

    memory->create(Pxx_fix,nmax,"fix:Pxx_fix");
    memory->create(Pxy_fix,nmax,"fix:Pxy_fix");
    memory->create(Pxz_fix,nmax,"fix:Pxz_fix");
    memory->create(Pyy_fix,nmax,"fix:Pyy_fix");
    memory->create(Pyz_fix,nmax,"fix:Pyz_fix");
    memory->create(Pzz_fix,nmax,"fix:Pzz_fix");

  }
  allocated = 1;
  // fprintf(screen,"\nCheckpoint 2");

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;

  // fprintf(screen,"\nCheckpoint 2B");

  // call widths and omega values from bespoke atom style.
  double *omega_SPH = atom->omega_SPH;
  double *width_SPH = atom->width_SPH;
  double *rho_SPH = atom->rho_SPH;
  double *u_SPH = atom->u_SPH;

  // fprintf(screen,"\nCheckpoint 2C");

  double SPH_weight_fact_self,SPH_weight_fact_ij,SPH_weight_fact_ji;

  double fx,fy,fz;

  double dt = update->dt;

  // fprintf(screen,"\nCheckpoint 2D");

  int *tag =atom->tag;

  double *mass = atom->mass;

  // fprintf(screen,"\nCheckpoint 2E");

  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int newton_pair = force->newton_pair;

  // fprintf(screen,"\nCheckpoint 2F");

  list = pair->list;

  // fprintf(screen,"\nCheckpoint 2G");
	
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  // fprintf(screen,"\nCheckpoint 2H");

  // // // fprintf(screen,"\nIn pair_bohm_dynamic_Mocz compute function...\n");

  // zero out per-atom arrays

  if (newton_pair) {
    for (i = 0; i < nall; i++){
      dxx_rho_fix[i] = 0.0;
      dxy_rho_fix[i] = 0.0;
      dxz_rho_fix[i] = 0.0;
      dyy_rho_fix[i] = 0.0;
      dyz_rho_fix[i] = 0.0;
      dzz_rho_fix[i] = 0.0;

      dx_rho_fix[i] = 0.0;
      dy_rho_fix[i] = 0.0;
      dz_rho_fix[i] = 0.0;

      Pxx_fix[i] = 0.0;
      Pxy_fix[i] = 0.0;
      Pxz_fix[i] = 0.0;
      Pyy_fix[i] = 0.0;
      Pyz_fix[i] = 0.0;
      Pzz_fix[i] = 0.0;
    }
  } 
  else{
    for (i = 0; i < nlocal; i++){
      dxx_rho_fix[i] = 0.0;
      dxy_rho_fix[i] = 0.0;
      dxz_rho_fix[i] = 0.0;
      dyy_rho_fix[i] = 0.0;
      dyz_rho_fix[i] = 0.0;
      dzz_rho_fix[i] = 0.0;

      dx_rho_fix[i] = 0.0;
      dy_rho_fix[i] = 0.0;
      dz_rho_fix[i] = 0.0;

      Pxx_fix[i] = 0.0;
      Pxy_fix[i] = 0.0;
      Pxz_fix[i] = 0.0;
      Pyy_fix[i] = 0.0;
      Pyz_fix[i] = 0.0;
      Pzz_fix[i] = 0.0;
    }
  }
  // fprintf(screen,"\nCheckpoint 3");

  // forward communicate velocities to ghost particles

  comm_forward = 3;
  comm_reverse = 3;
  commflag = 6;

  comm->forward_comm_fix(this);

  // fprintf(screen,"\nCheckpoint A");

  // loop over my atoms

  // 9 per-particle gradients to compute

  for (ii = 0; ii < inum; ii++) {

    // compute per-particle gradients for pressure tensor

    i = ilist[ii];

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    itype = type[i];
    jlist = firstneigh[i];

    jnum = numneigh[i];

    imass = mass[itype];

    h_i = width_SPH[i];
    h2_i = h_i*h_i;
    hm2_i = 1./h2_i;
    hm4_i = hm2_i*hm2_i;

    rho_i = rho_SPH[i];

    // 3D Gaussian prefactor
    gauss_pre_i = pi_fact*(1./(h_i*h_i*h_i));

    // no self interaction terms with Mocz type second order derivative

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

        delx_2 = delx*delx;
        dely_2 = dely*dely;
        delz_2 = delz*delz;

        m_gauss_ij = jmass*gauss_pre_i*exp(-(rsq)*hm2_i/2.);

        rho_j = rho_SPH[j];

        //apply Mocz factor to second order derivatives:

        dens_factor_ij = 1. - rho_i/rho_j;

        dx_rho_fix[i] += ((-delx)*hm2_i)*m_gauss_ij*dens_factor_ij;
        dy_rho_fix[i] += ((-dely)*hm2_i)*m_gauss_ij*dens_factor_ij;
        dz_rho_fix[i] += ((-delz)*hm2_i)*m_gauss_ij*dens_factor_ij;
        
        dxx_rho_fix[i] += hm2_i*(delx_2*hm2_i - 1.)*m_gauss_ij*dens_factor_ij;
        dxy_rho_fix[i] += hm4_i*(delx*dely)*m_gauss_ij*dens_factor_ij;
        dxz_rho_fix[i] += hm4_i*(delx*delz)*m_gauss_ij*dens_factor_ij;
        dyy_rho_fix[i] += hm2_i*(dely_2*hm2_i - 1.)*m_gauss_ij*dens_factor_ij;
        dyz_rho_fix[i] += hm4_i*(dely*delz)*m_gauss_ij*dens_factor_ij;
        dzz_rho_fix[i] += hm2_i*(delz_2*hm2_i - 1.)*m_gauss_ij*dens_factor_ij;

        if (newton_pair || j < nlocal) {

          h_j = width_SPH[j];
          h2_j = h_j*h_j;
          hm2_j = 1./h2_j;
          hm4_j = hm2_j*hm2_j;

          dens_factor_ji = 1. - rho_j/rho_i;

          gauss_pre_j = pi_fact*(1./(h_j*h_j*h_j));
        
          m_gauss_ji = imass*gauss_pre_j*exp(-(rsq)*hm2_j/2.);

          dx_rho_fix[j] += ((delx)*hm2_j)*m_gauss_ji*dens_factor_ji;
          dy_rho_fix[j] += ((dely)*hm2_j)*m_gauss_ji*dens_factor_ji;
          dz_rho_fix[j] += ((delz)*hm2_j)*m_gauss_ji*dens_factor_ji;

          dxx_rho_fix[j] += hm2_j*(delx_2*hm2_j - 1.)*m_gauss_ji*dens_factor_ji;
          dxy_rho_fix[j] += hm4_j*(delx*dely)*m_gauss_ji*dens_factor_ji;
          dxz_rho_fix[j] += hm4_j*(delx*delz)*m_gauss_ji*dens_factor_ji;
          dyy_rho_fix[j] += hm2_j*(dely_2*hm2_j - 1.)*m_gauss_ji*dens_factor_ji;
          dyz_rho_fix[j] += hm4_j*(dely*delz)*m_gauss_ji*dens_factor_ji;
          dzz_rho_fix[j] += hm2_j*(delz_2*hm2_j - 1.)*m_gauss_ji*dens_factor_ji;
        }
      }
    }
  }
  // fprintf(screen,"\nCheckpoint B");
  comm_forward = 9;
  comm_reverse = 9;
  commflag = 7;

  if (newton_pair) comm->reverse_comm_fix(this);

  comm->forward_comm_fix(this);

  // fprintf(screen,"\nCheckpoint C");

  for (ii = 0; ii < inum; ii++) {

    // compute per-particle gradients for pressure tensor

    i = ilist[ii];

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    itype = type[i];
    jlist = firstneigh[i];

    jnum = numneigh[i];

    imass = mass[itype];

    h_i = width_SPH[i];
    h2_i = h_i*h_i;
    hm2_i = 1./h2_i;

    // 3D Gaussian prefactor
    gauss_pre_i = pi_fact*(1./(h_i*h_i*h_i));

    // self contribution to Pressure Tensor

    SPH_weight_fact_self = (gauss_pre_i*imass)/rho_SPH[i];

    Pxx_fix[i] += SPH_weight_fact_self*gamma_factor*f_prefactor*((dx_rho_fix[i]*dx_rho_fix[i])/rho_SPH[i] - dxx_rho_fix[i]);
    Pxy_fix[i] += SPH_weight_fact_self*gamma_factor*f_prefactor*((dx_rho_fix[i]*dy_rho_fix[i])/rho_SPH[i] - dxy_rho_fix[i]);
    Pxz_fix[i] += SPH_weight_fact_self*gamma_factor*f_prefactor*((dx_rho_fix[i]*dz_rho_fix[i])/rho_SPH[i] - dxz_rho_fix[i]);
    Pyy_fix[i] += SPH_weight_fact_self*gamma_factor*f_prefactor*((dy_rho_fix[i]*dy_rho_fix[i])/rho_SPH[i] - dyy_rho_fix[i]);
    Pyz_fix[i] += SPH_weight_fact_self*gamma_factor*f_prefactor*((dy_rho_fix[i]*dz_rho_fix[i])/rho_SPH[i] - dyz_rho_fix[i]);
    Pzz_fix[i] += SPH_weight_fact_self*gamma_factor*f_prefactor*((dz_rho_fix[i]*dz_rho_fix[i])/rho_SPH[i] - dzz_rho_fix[i]);

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];

      rsq = delx*delx + dely*dely + delz*delz;

      jtype = type[j];
 
      if (rsq < cutsquared) {

        h_j = width_SPH[j];
        h2_j = h_j*h_j;
        hm2_j = 1./h2_j;

        jmass = mass[jtype];

        exp_ij = exp(-(rsq)*hm2_i/2.);

        gauss_pre_j = pi_fact*(1./(h_j*h_j*h_j));

        SPH_weight_fact_ij = (gauss_pre_i*jmass)/rho_SPH[j];

        Pxx_fix[i] += exp_ij*SPH_weight_fact_ij*gamma_factor*f_prefactor*((dx_rho_fix[j]*dx_rho_fix[j])/rho_SPH[j] - dxx_rho_fix[j]);
        Pxy_fix[i] += exp_ij*SPH_weight_fact_ij*gamma_factor*f_prefactor*((dx_rho_fix[j]*dy_rho_fix[j])/rho_SPH[j] - dxy_rho_fix[j]);
        Pxz_fix[i] += exp_ij*SPH_weight_fact_ij*gamma_factor*f_prefactor*((dx_rho_fix[j]*dz_rho_fix[j])/rho_SPH[j] - dxz_rho_fix[j]);
        Pyy_fix[i] += exp_ij*SPH_weight_fact_ij*gamma_factor*f_prefactor*((dy_rho_fix[j]*dy_rho_fix[j])/rho_SPH[j] - dyy_rho_fix[j]);
        Pyz_fix[i] += exp_ij*SPH_weight_fact_ij*gamma_factor*f_prefactor*((dy_rho_fix[j]*dz_rho_fix[j])/rho_SPH[j] - dyz_rho_fix[j]);
        Pzz_fix[i] += exp_ij*SPH_weight_fact_ij*gamma_factor*f_prefactor*((dz_rho_fix[j]*dz_rho_fix[j])/rho_SPH[j] - dzz_rho_fix[j]);

        if (newton_pair || j < nlocal) {

          exp_ji = exp(-(rsq)*hm2_j/2.);
          SPH_weight_fact_ji = (gauss_pre_j*imass)/rho_SPH[i];

          Pxx_fix[j] += exp_ji*SPH_weight_fact_ji*gamma_factor*f_prefactor*((dx_rho_fix[i]*dx_rho_fix[i])/rho_SPH[i] - dxx_rho_fix[i]);
          Pxy_fix[j] += exp_ji*SPH_weight_fact_ji*gamma_factor*f_prefactor*((dx_rho_fix[i]*dy_rho_fix[i])/rho_SPH[i] - dxy_rho_fix[i]);
          Pxz_fix[j] += exp_ji*SPH_weight_fact_ji*gamma_factor*f_prefactor*((dx_rho_fix[i]*dz_rho_fix[i])/rho_SPH[i] - dxz_rho_fix[i]);
          Pyy_fix[j] += exp_ji*SPH_weight_fact_ji*gamma_factor*f_prefactor*((dy_rho_fix[i]*dy_rho_fix[i])/rho_SPH[i] - dyy_rho_fix[i]);
          Pyz_fix[j] += exp_ji*SPH_weight_fact_ji*gamma_factor*f_prefactor*((dy_rho_fix[i]*dz_rho_fix[i])/rho_SPH[i] - dyz_rho_fix[i]);
          Pzz_fix[j] += exp_ji*SPH_weight_fact_ji*gamma_factor*f_prefactor*((dz_rho_fix[i]*dz_rho_fix[i])/rho_SPH[i] - dzz_rho_fix[i]);

        }
      }
    }
  }
  // fprintf(screen,"\nCheckpoint D");
  comm_forward = 6;
  comm_reverse = 6;
  commflag = 8;
  if (newton_pair) comm->reverse_comm_fix(this);
  comm->forward_comm_fix(this);
  // fprintf(screen,"\nCheckpoint E");

  for (ii = 0; ii < inum; ii++) {

    // compute per-particle gradients for pressure tensor

    i = ilist[ii];

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    vxtmp = v[i][0];
    vytmp = v[i][1];
    vztmp = v[i][2];

    itype = type[i];
    jlist = firstneigh[i];

    jnum = numneigh[i];

    imass = mass[itype];

    rho_i2 = rho_SPH[i]*rho_SPH[i];

    omega_i = omega_SPH[i];

    h_i = width_SPH[i];
    h2_i = h_i*h_i;
    hm2_i = 1./h2_i;

    // 3D Gaussian prefactor
    gauss_pre_i = pi_fact*(1./(h_i*h_i*h_i));

    u_prefact_i = (dt/(rho_i2*omega_i));

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];

      rsq = delx*delx + dely*dely + delz*delz;

      jtype = type[j];
 
      if (rsq < cutsquared) {

        delvx = vxtmp - v[j][0];
        delvy = vytmp - v[j][1];
        delvz = vztmp - v[j][2];

        h_j = width_SPH[j];
        h2_j = h_j*h_j;
        hm2_j = 1./h2_j;

        jmass = mass[jtype];
        ijmass = imass*jmass;

        exp_ij = exp(-(rsq)*hm2_i/2.);
        exp_ji = exp(-(rsq)*hm2_j/2.);

        omega_j = omega_SPH[j];

        gauss_pre_j = pi_fact*(1./(h_j*h_j*h_j));
      
        dx_Wij = gauss_pre_i*(-delx*hm2_i)*exp_ij;
        dy_Wij = gauss_pre_i*(-dely*hm2_i)*exp_ij;
        dz_Wij = gauss_pre_i*(-delz*hm2_i)*exp_ij;

        dx_Wji = gauss_pre_j*(-delx*hm2_j)*exp_ji;
        dy_Wji = gauss_pre_j*(-dely*hm2_j)*exp_ji;
        dz_Wji = gauss_pre_j*(-delz*hm2_j)*exp_ji;

        rho_j2 = rho_SPH[j]*rho_SPH[j];

        u_prefact_j = (dt/(rho_j2*omega_j));

        u_SPH[i] += u_prefact_i*ijmass*((Pxx_fix[i]*delvx + Pxy_fix[i]*delvy + Pxz_fix[i]*delvz)*dx_Wij + (Pxy_fix[i]*delvx + Pyy_fix[i]*delvy + Pyz_fix[i]*delvz)*dy_Wij + (Pxz_fix[i]*delvx + Pyz_fix[i]*delvy + Pzz_fix[i]*delvz)*dz_Wij);

        if (newton_pair || j < nlocal) {

          u_SPH[j] -= u_prefact_j*ijmass*((-Pxx_fix[j]*delvx - Pxy_fix[j]*delvy - Pxz_fix[j]*delvz)*dx_Wji + (-Pxy_fix[j]*delvx - Pyy_fix[j]*delvy - Pyz_fix[j]*delvz)*dy_Wji + (-Pxz_fix[j]*delvx - Pyz_fix[j]*delvy - Pzz_fix[j]*delvz)*dz_Wji);
        }
      }
    }
  }
  // fprintf(screen,"\nCheckpoint F");
  comm_forward = 1;
  comm_reverse = 1;
  commflag = 9;
  if (newton_pair) comm->reverse_comm_fix(this);
  comm->forward_comm_fix(this);
  // fprintf(screen,"\nCheckpoint G");

  // eflag_either = 1;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    eBohm += u_SPH[i];
    // if (eflag_atom) eatom[i] += u_SPH[i];
  }
  // fprintf(screen,"\nCheckpoint H");

}

int FixDynamicWidthsCoMBPMoczRobFtwo::pack_forward_comm(int n, int *list, double *buf,
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
  if (commflag == 6){
    double **v = atom->v;
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = v[j][0];
      buf[m++] = v[j][1];
      buf[m++] = v[j][2];
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
  if (commflag == 7){
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = dxx_rho_fix[j];
      buf[m++] = dxy_rho_fix[j];
      buf[m++] = dxz_rho_fix[j];
      buf[m++] = dyy_rho_fix[j];
      buf[m++] = dyz_rho_fix[j];
      buf[m++] = dzz_rho_fix[j];

      buf[m++] = dx_rho_fix[j];
      buf[m++] = dy_rho_fix[j];
      buf[m++] = dz_rho_fix[j];
    }
  }
  if (commflag == 8){
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = Pxx_fix[j];
      buf[m++] = Pxy_fix[j];
      buf[m++] = Pxz_fix[j];
      buf[m++] = Pyy_fix[j];
      buf[m++] = Pyz_fix[j];
      buf[m++] = Pzz_fix[j];
    }
  }
  if (commflag == 9){
    double *u_SPH = atom->u_SPH;
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = u_SPH[j];
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixDynamicWidthsCoMBPMoczRobFtwo::unpack_forward_comm(int n, int first, double *buf)
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
  if (commflag == 6){
    double **v = atom->v;
    for (i = first; i < last; i++){
      v[i][0] = buf[m++];
      v[i][1] = buf[m++];
      v[i][2] = buf[m++];
    }
  }
  if (commflag == 7){
    for (i = first; i < last; i++){
      dxx_rho_fix[i] = buf[m++];
      dxy_rho_fix[i] = buf[m++];
      dxz_rho_fix[i] = buf[m++];
      dyy_rho_fix[i] = buf[m++];
      dyz_rho_fix[i] = buf[m++];
      dzz_rho_fix[i] = buf[m++];

      dx_rho_fix[i] = buf[m++];
      dy_rho_fix[i] = buf[m++];
      dz_rho_fix[i] = buf[m++];
    }
  }
  if (commflag == 8){
    for (i = first; i < last; i++){
      Pxx_fix[i] = buf[m++];
      Pxy_fix[i] = buf[m++];
      Pxz_fix[i] = buf[m++];
      Pyy_fix[i] = buf[m++];
      Pyz_fix[i] = buf[m++];
      Pzz_fix[i] = buf[m++];
    }
  }
  if (commflag == 9){
    double *u_SPH = atom->u_SPH;
    for (i = first; i < last; i++){
      u_SPH[i] = buf[m++];
    }
  }
  
}

int FixDynamicWidthsCoMBPMoczRobFtwo::pack_reverse_comm(int n, int first, double *buf)
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
  if (commflag == 7){
    for (i = first; i < last; i++){
      buf[m++] = dxx_rho_fix[i];
      buf[m++] = dxy_rho_fix[i];
      buf[m++] = dxz_rho_fix[i];
      buf[m++] = dyy_rho_fix[i];
      buf[m++] = dyz_rho_fix[i];
      buf[m++] = dzz_rho_fix[i];

      buf[m++] = dx_rho_fix[i];
      buf[m++] = dy_rho_fix[i];
      buf[m++] = dz_rho_fix[i];
    }
  }
  if (commflag == 8){
    for (i = first; i < last; i++){
      buf[m++] = Pxx_fix[i];
      buf[m++] = Pxy_fix[i];
      buf[m++] = Pxz_fix[i];
      buf[m++] = Pyy_fix[i];
      buf[m++] = Pyz_fix[i];
      buf[m++] = Pzz_fix[i];
    }
  }
  if (commflag == 9){
    double *u_SPH = atom->u_SPH;
    // // // fprintf(screen,"In u_SPH pack reverse comm loop...");
    for (i = first; i < last; i++){
      buf[m++] = u_SPH[i];
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixDynamicWidthsCoMBPMoczRobFtwo::unpack_reverse_comm(int n, int *list, double *buf)
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
  if (commflag == 7){
    for (i = 0; i < n; i++) {
      j = list[i];
      dxx_rho_fix[j] += buf[m++];
      dxy_rho_fix[j] += buf[m++];
      dxz_rho_fix[j] += buf[m++];
      dyy_rho_fix[j] += buf[m++];
      dyz_rho_fix[j] += buf[m++];
      dzz_rho_fix[j] += buf[m++];

      dx_rho_fix[j] += buf[m++];
      dy_rho_fix[j] += buf[m++];
      dz_rho_fix[j] += buf[m++];
    }
  }
  if (commflag == 8){
    for (i = 0; i < n; i++) {
      j = list[i];
      Pxx_fix[j] += buf[m++];
      Pxy_fix[j] += buf[m++];
      Pxz_fix[j] += buf[m++];
      Pyy_fix[j] += buf[m++];
      Pyz_fix[j] += buf[m++];
      Pzz_fix[j] += buf[m++];
    }
  }
  if (commflag == 9){
    double *u_SPH = atom->u_SPH;
    // // // fprintf(screen,"In u_SPH unpack reverse comm loop...");
    for (i = 0; i < n; i++) {
      j = list[i];
      u_SPH[j] += buf[m++];
    }
  }
}

/* ----------------------------------------------------------------------
   allocate memory that depends on # of electrons
------------------------------------------------------------------------- */

void FixDynamicWidthsCoMBPMoczRobFtwo::allocate()
{
  // // fprintf(screen,"\nFixDynamicWidthsCoMBPMoczRobFtwo::allocate()\n");
  // Centre of Mass array x 0:N_ele, y N_ele:2*N_ele, z 2*N_ele:3_Nele
  mu = new double[3*N_electrons]; 
  mu_all = new double[3*N_electrons];
  tau = new double[3*N_electrons]; 
  tau_all = new double[3*N_electrons];
}

/* ----------------------------------------------------------------------
   deallocate memory that depends on # of electrons
------------------------------------------------------------------------- */

void FixDynamicWidthsCoMBPMoczRobFtwo::deallocate()
{
  // // fprintf(screen,"\nFixDynamicWidthsCoMBPMoczRobFtwo::deallocate()\n");
  delete [] mu;
  delete [] mu_all;
  delete [] tau;
  delete [] tau_all;
}

/* ----------------------------------------------------------------------
   potential energy SPH Bohmx
------------------------------------------------------------------------- */

double FixDynamicWidthsCoMBPMoczRobFtwo::compute_scalar()
{
  // only sum across procs one time

  if (eflag == 0) {
    MPI_Allreduce(&eBohm,&eBohm_all,1,MPI_DOUBLE,MPI_SUM,world);
    eflag = 1;
  }
  return eBohm_all;
}