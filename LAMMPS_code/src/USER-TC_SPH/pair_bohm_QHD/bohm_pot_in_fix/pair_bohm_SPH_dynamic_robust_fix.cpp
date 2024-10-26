/* -*- c++ -*- ----------------------------------------------------------
Bespoke pair_style to compute the Bohm force using an SPH-style pressure
tensor, with dynamic per-particle widths.

pair_bohm_SPH_dynamic_Mocz:
- Electrons as basic fluid element.
- Dynamic per-particle gaussian widths.
- Mocz 2015. form for second derivatives.
- Bohm potential calculated on total density of electron fluid (QHD Bohm Potential).

Thomas Campbell
------------------------------------------------------------------------- */
#include "pair_bohm_SPH_dynamic_robust_fix.h"
#include <mpi.h>
#include <cmath>
#include <cstring>
#include <math.h>
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"
#include "update.h"
#include "utils.h"

using namespace LAMMPS_NS;

#define MAXLINE 1024

/* ---------------------------------------------------------------------- */

PairBohmSPHDynamicRobustFix::PairBohmSPHDynamicRobustFix(LAMMPS *lmp) : Pair(lmp) {
  nmax = 0;

  manybody_flag = 1;

  dxx_rho = NULL;
  dxy_rho = NULL;
  dxz_rho = NULL;
  dyy_rho = NULL;
  dyz_rho = NULL;
  dzz_rho = NULL;

  // compute dx_rho terms internally in Bohm algorithm as F2 derivatives

  dx_rho = NULL;
  dy_rho = NULL;
  dz_rho = NULL;

  Pxx = NULL;
  Pxy = NULL;
  Pxz = NULL;
  Pyz = NULL;
  Pyy = NULL;
  Pzz = NULL;

  comm_forward = 9;
  comm_reverse = 9;

  hplanck  = force->hplanck;
  hbar = hplanck/(2*M_PI);
  // Bohm pressure prefactor
  f_prefactor = force->hhmrr2e * (hbar*hbar)/(4*e_mass);

  pi_fact = 1/pow(2*M_PI,1.5);
}

/* ---------------------------------------------------------------------- */

PairBohmSPHDynamicRobustFix::~PairBohmSPHDynamicRobustFix()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
    memory->destroy(dxx_rho);
    memory->destroy(dxy_rho);
    memory->destroy(dxz_rho);
    memory->destroy(dyy_rho);
    memory->destroy(dyz_rho);
    memory->destroy(dzz_rho);

    memory->destroy(dx_rho);
    memory->destroy(dy_rho);
    memory->destroy(dz_rho);

    memory->destroy(Pxx);
    memory->destroy(Pxy);
    memory->destroy(Pxz);
    memory->destroy(Pyy);
    memory->destroy(Pyz);
    memory->destroy(Pzz);
  }
}

/* ---------------------------------------------------------------------- */

void PairBohmSPHDynamicRobustFix::compute(int eflag, int vflag)
{
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
  double rho_i,rho_j;
  double bohm_pot;
  double u_prefact_i,u_prefact_j;

  // ev_init(eflag,vflag);

  int *ilist,*jlist,*numneigh,**firstneigh;

  if (atom->nmax > nmax) {
    // delete and create new memory arrays for any per-particle variables that need communicating.
    memory->destroy(dxx_rho);
    memory->destroy(dxy_rho);
    memory->destroy(dxz_rho);
    memory->destroy(dyy_rho);
    memory->destroy(dyz_rho);
    memory->destroy(dzz_rho);

    memory->destroy(dx_rho);
    memory->destroy(dy_rho);
    memory->destroy(dz_rho);

    memory->destroy(Pxx);
    memory->destroy(Pxy);
    memory->destroy(Pxz);
    memory->destroy(Pyy);
    memory->destroy(Pyz);
    memory->destroy(Pzz);
    
    nmax = atom->nmax;

    memory->create(dxx_rho,nmax,"pair:dxx_rho");
    memory->create(dxy_rho,nmax,"pair:dxy_rho");
    memory->create(dxz_rho,nmax,"pair:dxz_rho");
    memory->create(dyy_rho,nmax,"pair:dyy_rho");
    memory->create(dyz_rho,nmax,"pair:dyz_rho");
    memory->create(dzz_rho,nmax,"pair:dzz_rho");

    memory->create(dx_rho,nmax,"pair:dx_rho");
    memory->create(dy_rho,nmax,"pair:dy_rho");
    memory->create(dz_rho,nmax,"pair:dz_rho");

    memory->create(Pxx,nmax,"pair:Pxx");
    memory->create(Pxy,nmax,"pair:Pxy");
    memory->create(Pxz,nmax,"pair:Pxz");
    memory->create(Pyy,nmax,"pair:Pyy");
    memory->create(Pyz,nmax,"pair:Pyz");
    memory->create(Pzz,nmax,"pair:Pzz");

  }

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;

  // call widths and omega values from bespoke atom style.
  double *omega_SPH = atom->omega_SPH;
  double *width_SPH = atom->width_SPH;
  double *rho_SPH = atom->rho_SPH;
  double *u_SPH = atom->u_SPH;

  double SPH_weight_fact_self,SPH_weight_fact_ij,SPH_weight_fact_ji;

  double fx,fy,fz;

  double dt = update->dt;

  int *tag =atom->tag;

  double *mass = atom->mass;

  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // // fprintf(screen,"\nIn pair_bohm_dynamic_Mocz compute function...\n");

  // zero out per-atom arrays

  if (newton_pair) {
    for (i = 0; i < nall; i++){
      dxx_rho[i] = 0.0;
      dxy_rho[i] = 0.0;
      dxz_rho[i] = 0.0;
      dyy_rho[i] = 0.0;
      dyz_rho[i] = 0.0;
      dzz_rho[i] = 0.0;

      dx_rho[i] = 0.0;
      dy_rho[i] = 0.0;
      dz_rho[i] = 0.0;

      Pxx[i] = 0.0;
      Pxy[i] = 0.0;
      Pxz[i] = 0.0;
      Pyy[i] = 0.0;
      Pyz[i] = 0.0;
      Pzz[i] = 0.0;
    }
  } 
  else{
    for (i = 0; i < nlocal; i++){
      dxx_rho[i] = 0.0;
      dxy_rho[i] = 0.0;
      dxz_rho[i] = 0.0;
      dyy_rho[i] = 0.0;
      dyz_rho[i] = 0.0;
      dzz_rho[i] = 0.0;

      dx_rho[i] = 0.0;
      dy_rho[i] = 0.0;
      dz_rho[i] = 0.0;

      Pxx[i] = 0.0;
      Pxy[i] = 0.0;
      Pxz[i] = 0.0;
      Pyy[i] = 0.0;
      Pyz[i] = 0.0;
      Pzz[i] = 0.0;
    }
  }

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

        dx_rho[i] += ((-delx)*hm2_i)*m_gauss_ij;
        dy_rho[i] += ((-dely)*hm2_i)*m_gauss_ij;
        dz_rho[i] += ((-delz)*hm2_i)*m_gauss_ij;
        
        dxx_rho[i] += hm2_i*(delx_2*hm2_i - 1.)*m_gauss_ij;
        dxy_rho[i] += hm4_i*(delx*dely)*m_gauss_ij;
        dxz_rho[i] += hm4_i*(delx*delz)*m_gauss_ij;
        dyy_rho[i] += hm2_i*(dely_2*hm2_i - 1.)*m_gauss_ij;
        dyz_rho[i] += hm4_i*(dely*delz)*m_gauss_ij;
        dzz_rho[i] += hm2_i*(delz_2*hm2_i - 1.)*m_gauss_ij;

        if (newton_pair || j < nlocal) {

          h_j = width_SPH[j];
          h2_j = h_j*h_j;
          hm2_j = 1./h2_j;
          hm4_j = hm2_j*hm2_j;

          gauss_pre_j = pi_fact*(1./(h_j*h_j*h_j));
        
          m_gauss_ji = imass*gauss_pre_j*exp(-(rsq)*hm2_j/2.);

          dx_rho[j] += ((delx)*hm2_j)*m_gauss_ji;
          dy_rho[j] += ((dely)*hm2_j)*m_gauss_ji;
          dz_rho[j] += ((delz)*hm2_j)*m_gauss_ji;

          dxx_rho[j] += hm2_j*(delx_2*hm2_j - 1.)*m_gauss_ji;
          dxy_rho[j] += hm4_j*(delx*dely)*m_gauss_ji;
          dxz_rho[j] += hm4_j*(delx*delz)*m_gauss_ji;
          dyy_rho[j] += hm2_j*(dely_2*hm2_j - 1.)*m_gauss_ji;
          dyz_rho[j] += hm4_j*(dely*delz)*m_gauss_ji;
          dzz_rho[j] += hm2_j*(delz_2*hm2_j - 1.)*m_gauss_ji;
        }
      }
    }
  }
  comm_forward = 9;
  comm_reverse = 9;
  commflag = 0;

  if (newton_pair) comm->reverse_comm_pair(this);

  comm->forward_comm_pair(this);

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

    Pxx[i] += SPH_weight_fact_self*gamma_factor*f_prefactor*((dx_rho[i]*dx_rho[i])/rho_SPH[i] - dxx_rho[i]);
    Pxy[i] += SPH_weight_fact_self*gamma_factor*f_prefactor*((dx_rho[i]*dy_rho[i])/rho_SPH[i] - dxy_rho[i]);
    Pxz[i] += SPH_weight_fact_self*gamma_factor*f_prefactor*((dx_rho[i]*dz_rho[i])/rho_SPH[i] - dxz_rho[i]);
    Pyy[i] += SPH_weight_fact_self*gamma_factor*f_prefactor*((dy_rho[i]*dy_rho[i])/rho_SPH[i] - dyy_rho[i]);
    Pyz[i] += SPH_weight_fact_self*gamma_factor*f_prefactor*((dy_rho[i]*dz_rho[i])/rho_SPH[i] - dyz_rho[i]);
    Pzz[i] += SPH_weight_fact_self*gamma_factor*f_prefactor*((dz_rho[i]*dz_rho[i])/rho_SPH[i] - dzz_rho[i]);

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

        Pxx[i] += exp_ij*SPH_weight_fact_ij*gamma_factor*f_prefactor*((dx_rho[j]*dx_rho[j])/rho_SPH[j] - dxx_rho[j]);
        Pxy[i] += exp_ij*SPH_weight_fact_ij*gamma_factor*f_prefactor*((dx_rho[j]*dy_rho[j])/rho_SPH[j] - dxy_rho[j]);
        Pxz[i] += exp_ij*SPH_weight_fact_ij*gamma_factor*f_prefactor*((dx_rho[j]*dz_rho[j])/rho_SPH[j] - dxz_rho[j]);
        Pyy[i] += exp_ij*SPH_weight_fact_ij*gamma_factor*f_prefactor*((dy_rho[j]*dy_rho[j])/rho_SPH[j] - dyy_rho[j]);
        Pyz[i] += exp_ij*SPH_weight_fact_ij*gamma_factor*f_prefactor*((dy_rho[j]*dz_rho[j])/rho_SPH[j] - dyz_rho[j]);
        Pzz[i] += exp_ij*SPH_weight_fact_ij*gamma_factor*f_prefactor*((dz_rho[j]*dz_rho[j])/rho_SPH[j] - dzz_rho[j]);

        if (newton_pair || j < nlocal) {

          exp_ji = exp(-(rsq)*hm2_j/2.);
          SPH_weight_fact_ji = (gauss_pre_j*imass)/rho_SPH[i];

          Pxx[j] += exp_ji*SPH_weight_fact_ji*gamma_factor*f_prefactor*((dx_rho[i]*dx_rho[i])/rho_SPH[i] - dxx_rho[i]);
          Pxy[j] += exp_ji*SPH_weight_fact_ji*gamma_factor*f_prefactor*((dx_rho[i]*dy_rho[i])/rho_SPH[i] - dxy_rho[i]);
          Pxz[j] += exp_ji*SPH_weight_fact_ji*gamma_factor*f_prefactor*((dx_rho[i]*dz_rho[i])/rho_SPH[i] - dxz_rho[i]);
          Pyy[j] += exp_ji*SPH_weight_fact_ji*gamma_factor*f_prefactor*((dy_rho[i]*dy_rho[i])/rho_SPH[i] - dyy_rho[i]);
          Pyz[j] += exp_ji*SPH_weight_fact_ji*gamma_factor*f_prefactor*((dy_rho[i]*dz_rho[i])/rho_SPH[i] - dyz_rho[i]);
          Pzz[j] += exp_ji*SPH_weight_fact_ji*gamma_factor*f_prefactor*((dz_rho[i]*dz_rho[i])/rho_SPH[i] - dzz_rho[i]);

        }
      }
    }
  }
  comm_forward = 6;
  comm_reverse = 6;
  commflag = 2;
  if (newton_pair) comm->reverse_comm_pair(this);
  comm->forward_comm_pair(this);

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
        // compute force terms from pressure tensor here

        fx = -ijmass*((Pxx[i]*dx_Wij + Pxy[i]*dy_Wij + Pxz[i]*dz_Wij)/(rho_i2*omega_i) + (Pxx[j]*dx_Wji + Pxy[j]*dy_Wji + Pxz[j]*dz_Wji)/(rho_j2*omega_j));
        fy = -ijmass*((Pxy[i]*dx_Wij + Pyy[i]*dy_Wij + Pyz[i]*dz_Wij)/(rho_i2*omega_i) + (Pxy[j]*dx_Wji + Pyy[j]*dy_Wji + Pyz[j]*dz_Wji)/(rho_j2*omega_j));
        fz = -ijmass*((Pxz[i]*dx_Wij + Pyz[i]*dy_Wij + Pzz[i]*dz_Wij)/(rho_i2*omega_i) + (Pxz[j]*dx_Wji + Pyz[j]*dy_Wji + Pzz[j]*dz_Wji)/(rho_j2*omega_j));
        
        // add flags to check if force terms become attractive 

        // if (delx/fx < 0) {
        //   fprintf(screen,"\n Particle ID %d attractive Bohm x force",tag[i]);
        //   fprintf(screen,"\n i fx = %16.16f",fx);
        //   fprintf(screen,"\n xi - xj = %16.16f",delx);
        // }
        // if (dely/fy < 0) {
        //   fprintf(screen,"\n Particle ID %d attractive Bohm y force",tag[i]);
        //   fprintf(screen,"\n i fy = %16.16f",fy);
        //   fprintf(screen,"\n yi - yj = %16.16f",dely);
        // }
        // if (delz/fz < 0) {
        //   fprintf(screen,"\n Particle ID %d attractive Bohm z force",tag[i]);
        //   fprintf(screen,"\n i fz = %16.16f",fz);
        //   fprintf(screen,"\n zi - zj = %16.16f",delz);
        // }

        // if (delx/fx < 0) {
        //   if (dely/fy < 0) {
        //     if (delz/fz < 0) {
        //       fprintf(screen,"\n Particle IDs i %d j %d attractive Bohm force",tag[i],tag[j]);
        //       fprintf(screen,"\n i fx = %16.16f",fx);
        //       fprintf(screen,"\n xi - xj = %16.16f",delx);
        //       fprintf(screen,"\n i fy = %16.16f",fy);
        //       fprintf(screen,"\n yi - yj = %16.16f",dely);
        //       fprintf(screen,"\n i fz = %16.16f",fz);
        //       fprintf(screen,"\n zi - zj = %16.16f",delz);
        //     }
        //   }
        // }

        f[i][0] += fx;
        f[i][1] += fy;
        f[i][2] += fz;

        // fprintf(screen,"\nBohm fx = %16.16f",fx);

        if (newton_pair || j < nlocal) {
          f[j][0] -= fx;
          f[j][1] -= fy;
          f[j][2] -= fz;

        }

        eflag_either = 0;
        if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,0.0,0.0,
                        fx,fy,fz,delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairBohmSPHDynamicRobustFix::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut,n+1,n+1,"pair:cut");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairBohmSPHDynamicRobustFix::settings(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR,"Illegal pair_style command. Require 2 input arguments.");

  if (comm->ghost_velocity != 1) error->all(FLERR,"Illegal pair_style command. Require ghost atom velocity.");

  cut_global = force->numeric(FLERR,arg[0]);
  gamma_factor = force->numeric(FLERR,arg[1]);

  cutsquared = cut_global*cut_global;

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairBohmSPHDynamicRobustFix::coeff(int narg, char **arg)
{
  if (narg < 2 || narg > 3)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

// Specify a different cutoff in pair_coeff step:

  double cut_one = cut_global;
  if (narg == 3) cut_one = force->numeric(FLERR,arg[2]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}


/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairBohmSPHDynamicRobustFix::init_style()
{
  neighbor->request(this,instance_me);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairBohmSPHDynamicRobustFix::init_one(int i, int j)
{
  if (setflag[i][j] == 0)
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);

  return cut[i][j];
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairBohmSPHDynamicRobustFix::write_restart(FILE *fp)
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

void PairBohmSPHDynamicRobustFix::read_restart(FILE *fp)
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

void PairBohmSPHDynamicRobustFix::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairBohmSPHDynamicRobustFix::read_restart_settings(FILE *fp)
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

int PairBohmSPHDynamicRobustFix::pack_forward_comm(int n, int *list, double *buf,
                               int /*pbc_flag*/, int * /*pbc*/)
{
  int i,j,m;

  m = 0;
  if (commflag == 0){
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = dxx_rho[j];
      buf[m++] = dxy_rho[j];
      buf[m++] = dxz_rho[j];
      buf[m++] = dyy_rho[j];
      buf[m++] = dyz_rho[j];
      buf[m++] = dzz_rho[j];

      buf[m++] = dx_rho[j];
      buf[m++] = dy_rho[j];
      buf[m++] = dz_rho[j];
    }
  }
  if (commflag == 1){
    double *u_SPH = atom->u_SPH;
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = u_SPH[j];
    }
  }
  if (commflag == 2){
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = Pxx[j];
      buf[m++] = Pxy[j];
      buf[m++] = Pxz[j];
      buf[m++] = Pyy[j];
      buf[m++] = Pyz[j];
      buf[m++] = Pzz[j];
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairBohmSPHDynamicRobustFix::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  if (commflag == 0){
    for (i = first; i < last; i++){
      dxx_rho[i] = buf[m++];
      dxy_rho[i] = buf[m++];
      dxz_rho[i] = buf[m++];
      dyy_rho[i] = buf[m++];
      dyz_rho[i] = buf[m++];
      dzz_rho[i] = buf[m++];

      dx_rho[i] = buf[m++];
      dy_rho[i] = buf[m++];
      dz_rho[i] = buf[m++];
    }
  }
  if (commflag == 1){
    double *u_SPH = atom->u_SPH;
    for (i = first; i < last; i++){
      u_SPH[i] = buf[m++];
    }
  }
  if (commflag == 2){
    for (i = first; i < last; i++){
      Pxx[i] = buf[m++];
      Pxy[i] = buf[m++];
      Pxz[i] = buf[m++];
      Pyy[i] = buf[m++];
      Pyz[i] = buf[m++];
      Pzz[i] = buf[m++];
    }
  }
  
}

/* ---------------------------------------------------------------------- */

int PairBohmSPHDynamicRobustFix::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  if (commflag == 0){
    for (i = first; i < last; i++){
      buf[m++] = dxx_rho[i];
      buf[m++] = dxy_rho[i];
      buf[m++] = dxz_rho[i];
      buf[m++] = dyy_rho[i];
      buf[m++] = dyz_rho[i];
      buf[m++] = dzz_rho[i];

      buf[m++] = dx_rho[i];
      buf[m++] = dy_rho[i];
      buf[m++] = dz_rho[i];
    }
  }
  if (commflag == 1){
    double *u_SPH = atom->u_SPH;
    // // fprintf(screen,"In u_SPH pack reverse comm loop...");
    for (i = first; i < last; i++){
      buf[m++] = u_SPH[i];
    }
  }
  if (commflag == 2){
    for (i = first; i < last; i++){
      buf[m++] = Pxx[i];
      buf[m++] = Pxy[i];
      buf[m++] = Pxz[i];
      buf[m++] = Pyy[i];
      buf[m++] = Pyz[i];
      buf[m++] = Pzz[i];
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairBohmSPHDynamicRobustFix::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  if (commflag == 0){
    for (i = 0; i < n; i++) {
      j = list[i];
      dxx_rho[j] += buf[m++];
      dxy_rho[j] += buf[m++];
      dxz_rho[j] += buf[m++];
      dyy_rho[j] += buf[m++];
      dyz_rho[j] += buf[m++];
      dzz_rho[j] += buf[m++];

      dx_rho[j] += buf[m++];
      dy_rho[j] += buf[m++];
      dz_rho[j] += buf[m++];
    }
  }
  if (commflag == 1){
    double *u_SPH = atom->u_SPH;
    // // fprintf(screen,"In u_SPH unpack reverse comm loop...");
    for (i = 0; i < n; i++) {
      j = list[i];
      u_SPH[j] += buf[m++];
    }
  }
  if (commflag == 2){
    for (i = 0; i < n; i++) {
      j = list[i];
      Pxx[j] += buf[m++];
      Pxy[j] += buf[m++];
      Pxz[j] += buf[m++];
      Pyy[j] += buf[m++];
      Pyz[j] += buf[m++];
      Pzz[j] += buf[m++];
    }
  }
}

/* ---------------------------------------------------------------------- */

void *PairBohmSPHDynamicRobustFix::extract(const char *str, int &dim)
{
  dim = 2;
  return NULL;
}
