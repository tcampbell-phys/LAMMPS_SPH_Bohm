/* -*- c++ -*- ----------------------------------------------------------
Bespoke pair_style to compute the Bohm force using an SPH-style pressure
tensor, with dynamic per-particle widths.

pair_bohm_SPH_dynamic_Mocz:
- Electrons as basic fluid element.
- Dynamic per-particle gaussian widths.
- Mocz 2015. form for second derivatives.

Thomas Campbell
------------------------------------------------------------------------- */
#include "pair_bohm_SPH_dynamic_Mocz.h"
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

PairBohmSPHDynamicMocz::PairBohmSPHDynamicMocz(LAMMPS *lmp) : Pair(lmp) {
  nmax = 0;

  manybody_flag = 1;

  dx_rho = NULL;
  dy_rho = NULL;
  dz_rho = NULL;
  dxx_rho = NULL;
  dxy_rho = NULL;
  dxz_rho = NULL;
  dyy_rho = NULL;
  dyz_rho = NULL;
  dzz_rho = NULL;

  comm_forward = 9;
  comm_reverse = 9;
}

/* ---------------------------------------------------------------------- */

PairBohmSPHDynamicMocz::~PairBohmSPHDynamicMocz()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
    memory->destroy(dx_rho);
    memory->destroy(dy_rho);
    memory->destroy(dz_rho);
    memory->destroy(dxx_rho);
    memory->destroy(dxy_rho);
    memory->destroy(dxz_rho);
    memory->destroy(dyy_rho);
    memory->destroy(dyz_rho);
    memory->destroy(dzz_rho);
  }
}

/* ---------------------------------------------------------------------- */

void PairBohmSPHDynamicMocz::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz;
  double xtmp2,ytmp2,ztmp2,delx2,dely2,delz2;
  double delx_2,dely_2,delz_2,rsq;
  double pi_fact;
  double gauss_pre_i;
  double gauss_pre_j;
  double h_i,h2_i,hm2_i,hm4_i;
  double h_j,h2_j,hm2_j,hm4_j;
  double imass,jmass,ijmass;
  double m_gauss_ij,m_gauss_ji;
  double exp_ij,exp_ji;
  double omega_i,omega_j;
  double Pixx,Pixy,Pixz,Piyy,Piyz,Pizz;
  double Pjxx,Pjxy,Pjxz,Pjyy,Pjyz,Pjzz;
  double dx_Wij,dy_Wij,dz_Wij;
  double dx_Wji,dy_Wji,dz_Wji;
  double rho_i2,rho_j2;
  double dens_factor_ij,dens_factor_ji;
  double rho_i,rho_j;
  double bohm_pot;
  double cutsquared;
  
  hplanck  = force->hplanck;
  hbar = hplanck/(2*M_PI);

  // Bohm pressure prefactor
  f_prefactor = (hbar*hbar)/(4*e_mass);

  int *ilist,*jlist,*numneigh,**firstneigh;

  ev_init(eflag,vflag);

  if (atom->nmax > nmax) {
    // delete and create new memory arrays for any per-particle variables that need communicating.
    memory->destroy(dx_rho);
    memory->destroy(dy_rho);
    memory->destroy(dz_rho);
    memory->destroy(dxx_rho);
    memory->destroy(dxy_rho);
    memory->destroy(dxz_rho);
    memory->destroy(dyy_rho);
    memory->destroy(dyz_rho);
    memory->destroy(dzz_rho);
    
    nmax = atom->nmax;

    memory->create(dx_rho,nmax,"pair:dx_rho");
    memory->create(dy_rho,nmax,"pair:dy_rho");
    memory->create(dz_rho,nmax,"pair:dz_rho");
    memory->create(dxx_rho,nmax,"pair:dxx_rho");
    memory->create(dxy_rho,nmax,"pair:dxy_rho");
    memory->create(dxz_rho,nmax,"pair:dxz_rho");
    memory->create(dyy_rho,nmax,"pair:dyy_rho");
    memory->create(dyz_rho,nmax,"pair:dyz_rho");
    memory->create(dzz_rho,nmax,"pair:dzz_rho");

  }

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;

  // call widths and omega values from bespoke atom style.
  double *omega_SPH = atom->omega_SPH;
  double *width_SPH = atom->width_SPH;
  double *rho_SPH = atom->rho_SPH;

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

  cutsquared = cut_global*cut_global;

  // zero out per-atom arrays

  if (newton_pair) {
    for (i = 0; i < nall; i++){
      dx_rho[i] = 0.0;
      dy_rho[i] = 0.0;
      dz_rho[i] = 0.0;
      dxx_rho[i] = 0.0;
      dxy_rho[i] = 0.0;
      dxz_rho[i] = 0.0;
      dyy_rho[i] = 0.0;
      dyz_rho[i] = 0.0;
      dzz_rho[i] = 0.0;
    }
  } 
  else{
    for (i = 0; i < nlocal; i++){
      dx_rho[i] = 0.0;
      dy_rho[i] = 0.0;
      dz_rho[i] = 0.0;
      dxx_rho[i] = 0.0;
      dxy_rho[i] = 0.0;
      dxz_rho[i] = 0.0;
      dyy_rho[i] = 0.0;
      dyz_rho[i] = 0.0;
      dzz_rho[i] = 0.0;
    }
  }

  // loop over my atoms

  // 9 per-particle gradients to compute

  pi_fact = 1/pow(2*M_PI,1.5);

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

    // fprintf(screen,"break_statement A\n");
    // fprintf(screen,"Particle at...\n");
    // fprintf(screen,"x[i] = %.12f \n",xtmp);
    // fprintf(screen,"y[i] = %.12f \n",ytmp);
    // fprintf(screen,"z[i] = %.12f \n",ztmp);

    h_i = width_SPH[i];
    h2_i = h_i*h_i;
    hm2_i = 1./h2_i;
    hm4_i = hm2_i*hm2_i;

    rho_i = rho_SPH[i];

    // 3D Gaussian prefactor
    gauss_pre_i = pi_fact*(1./(h_i*h_i*h_i));

    // no self interaction terms with Mocz type second order derivative

    for (jj = 0; jj < jnum; jj++) {
      // fprintf(screen,"break_statement B\n");
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];

      rsq = delx*delx + dely*dely + delz*delz;

      jtype = type[j];
 
      if (rsq < cutsquared) {
        // fprintf(screen,"break_statement C\n");

        // fprintf(screen,"Neighbour at...\n");
        // fprintf(screen,"x[j] = %.12f \n",x[j][0]);
        // fprintf(screen,"y[j] = %.12f \n",x[j][1]);
        // fprintf(screen,"z[j] = %.12f \n",x[j][2]);

        jmass = mass[jtype];

        delx_2 = delx*delx;
        dely_2 = dely*dely;
        delz_2 = delz*delz;

        m_gauss_ij = jmass*gauss_pre_i*exp(-(rsq)*hm2_i/2.);

        rho_j = rho_SPH[j];

        dens_factor_ij = 1. - rho_i/rho_j;
        
        dx_rho[i] += ((-delx)*hm2_i)*m_gauss_ij;
        dy_rho[i] += ((-dely)*hm2_i)*m_gauss_ij;
        dz_rho[i] += ((-delz)*hm2_i)*m_gauss_ij;
        dxx_rho[i] += hm2_i*(delx_2*hm2_i - 1.)*m_gauss_ij*dens_factor_ij;
        dxy_rho[i] += hm4_i*(delx*dely)*m_gauss_ij*dens_factor_ij;
        dxz_rho[i] += hm4_i*(delx*delz)*m_gauss_ij*dens_factor_ij;
        dyy_rho[i] += hm2_i*(dely_2*hm2_i - 1.)*m_gauss_ij*dens_factor_ij;
        dyz_rho[i] += hm4_i*(dely*delz)*m_gauss_ij*dens_factor_ij;
        dzz_rho[i] += hm2_i*(delz_2*hm2_i - 1.)*m_gauss_ij*dens_factor_ij;

        if (newton_pair || j < nlocal) {
          // fprintf(screen,"break_statement D\n");

          h_j = width_SPH[j];
          h2_j = h_j*h_j;
          hm2_j = 1./h2_j;
          hm4_j = hm2_j*hm2_j;

          dens_factor_ji = 1. - rho_j/rho_i;

          gauss_pre_j = pi_fact*(1./(h_j*h_j*h_j));
        
          m_gauss_ji = imass*gauss_pre_j*exp(-(rsq)*hm2_j/2.);

          dx_rho[j] += ((delx)*hm2_j)*m_gauss_ji;
          dy_rho[j] += ((dely)*hm2_j)*m_gauss_ji;
          dz_rho[j] += ((delz)*hm2_j)*m_gauss_ji;
          dxx_rho[j] += hm2_j*(delx_2*hm2_j - 1.)*m_gauss_ji*dens_factor_ji;
          dxy_rho[j] += hm4_j*(delx*dely)*m_gauss_ji*dens_factor_ji;
          dxz_rho[j] += hm4_j*(delx*delz)*m_gauss_ji*dens_factor_ji;
          dyy_rho[j] += hm2_j*(dely_2*hm2_j - 1.)*m_gauss_ji*dens_factor_ji;
          dyz_rho[j] += hm4_j*(dely*delz)*m_gauss_ji*dens_factor_ji;
          dzz_rho[j] += hm2_j*(delz_2*hm2_j - 1.)*m_gauss_ji*dens_factor_ji;
        }
      }
    }
  }
  
  commflag = 0;

  if (newton_pair) comm->reverse_comm_pair(this);

  comm->forward_comm_pair(this);

  for (ii = 0; ii < inum; ii++) {
    // fprintf(screen,"break_statement E\n");

    // compute per-particle gradients for pressure tensor

    i = ilist[ii];

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    // fprintf(screen,"\nParticle at...\n");
    // fprintf(screen,"x =  %.12f\n", xtmp);
    // fprintf(screen,"y =  %.12f\n", ytmp);
    // fprintf(screen,"z =  %.12f\n", ztmp);

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

    Pixx = gamma_factor*f_prefactor*((dx_rho[i]*dx_rho[i])/rho_SPH[i] - dxx_rho[i]);
    Pixy = gamma_factor*f_prefactor*((dx_rho[i]*dy_rho[i])/rho_SPH[i] - dxy_rho[i]);
    Pixz = gamma_factor*f_prefactor*((dx_rho[i]*dz_rho[i])/rho_SPH[i] - dxz_rho[i]);
    Piyy = gamma_factor*f_prefactor*((dy_rho[i]*dy_rho[i])/rho_SPH[i] - dyy_rho[i]);
    Piyz = gamma_factor*f_prefactor*((dy_rho[i]*dz_rho[i])/rho_SPH[i] - dyz_rho[i]);
    Pizz = gamma_factor*f_prefactor*((dz_rho[i]*dz_rho[i])/rho_SPH[i] - dzz_rho[i]);

    // fprintf(screen,"Pixx =  %.12f\n", Pixx);
    // fprintf(screen,"Pixy =  %.12f\n", Pixy);
    // fprintf(screen,"Pixz =  %.12f\n", Pixz);
    // fprintf(screen,"Piyy =  %.12f\n", Piyy);
    // fprintf(screen,"Piyz =  %.12f\n", Piyz);
    // fprintf(screen,"Pizz =  %.12f\n", Pizz);
    // fprintf(screen,"width[i] =  %.12f\n", width_SPH[i]);
    // fprintf(screen,"omega_i =  %.12f\n", omega_i);
    // fprintf(screen,"rho[i] =  %.12f\n", rho_SPH[i]);
    

    // Bohm potential calculation
    bohm_pot = -gamma_factor*f_prefactor*((dxx_rho[i] + dyy_rho[i]+ dzz_rho[i])/rho_SPH[i] - (dx_rho[i]*dx_rho[i] + dy_rho[i]*dy_rho[i] + dz_rho[i]*dz_rho[i])/(2*rho_i2));
    if (eflag_global) eng_vdwl += bohm_pot;

    for (jj = 0; jj < jnum; jj++) {
      // fprintf(screen,"break_statement F\n");
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];

      rsq = delx*delx + dely*dely + delz*delz;

      jtype = type[j];
 
      if (rsq < cutsquared) {
        // fprintf(screen,"break_statement G\n");

        // fprintf(screen,"neighbour at...\n");
        // fprintf(screen,"x =  %.12f\n", x[j][0]);
        // fprintf(screen,"y =  %.12f\n", x[j][1]);
        // fprintf(screen,"z =  %.12f\n", x[j][2]);

        h_j = width_SPH[j];
        h2_j = h_j*h_j;
        hm2_j = 1./h2_j;

        exp_ij = exp(-(rsq)*hm2_i/2.);
        exp_ji = exp(-(rsq)*hm2_j/2.);

        omega_j = omega_SPH[j];

        gauss_pre_j = pi_fact*(1./(h_j*h_j*h_j));

        Pjxx = gamma_factor*f_prefactor*((dx_rho[j]*dx_rho[j])/rho_SPH[j] - dxx_rho[j]);
        Pjxy = gamma_factor*f_prefactor*((dx_rho[j]*dy_rho[j])/rho_SPH[j] - dxy_rho[j]);
        Pjxz = gamma_factor*f_prefactor*((dx_rho[j]*dz_rho[j])/rho_SPH[j] - dxz_rho[j]);
        Pjyy = gamma_factor*f_prefactor*((dy_rho[j]*dy_rho[j])/rho_SPH[j] - dyy_rho[j]);
        Pjyz = gamma_factor*f_prefactor*((dy_rho[j]*dz_rho[j])/rho_SPH[j] - dyz_rho[j]);
        Pjzz = gamma_factor*f_prefactor*((dz_rho[j]*dz_rho[j])/rho_SPH[j] - dzz_rho[j]);
      
        dx_Wij = gauss_pre_i*(-delx*hm2_i)*exp_ij;
        dy_Wij = gauss_pre_i*(-dely*hm2_i)*exp_ij;
        dz_Wij = gauss_pre_i*(-delz*hm2_i)*exp_ij;

        dx_Wji = gauss_pre_j*(-delx*hm2_j)*exp_ji;
        dy_Wji = gauss_pre_j*(-dely*hm2_j)*exp_ji;
        dz_Wji = gauss_pre_j*(-delz*hm2_j)*exp_ji;

        // fprintf(screen,"Pjxx =  %.12f\n", Pjxx);
        // fprintf(screen,"Pjxy =  %.12f\n", Pjxy);
        // fprintf(screen,"Pjxz =  %.12f\n", Pjxz);
        // fprintf(screen,"Pjyy =  %.12f\n", Pjyy);
        // fprintf(screen,"Pjyz =  %.12f\n", Pjyz);
        // fprintf(screen,"Pjzz =  %.12f\n", Pjzz);
        // fprintf(screen,"width[j] =  %.12f\n", width_SPH[j]);
        // fprintf(screen,"omega_j =  %.12f\n", omega_j);
        // fprintf(screen,"rho[j] =  %.12f\n", rho_SPH[j]);
        // fprintf(screen,"dx_Wij =  %.12f\n", dx_Wij);
        // fprintf(screen,"dy_Wij =  %.12f\n", dy_Wij);
        // fprintf(screen,"dz_Wij =  %.12f\n", dz_Wij);
        // fprintf(screen,"dx_Wji =  %.12f\n", dx_Wji);
        // fprintf(screen,"dy_Wji =  %.12f\n", dy_Wji);
        // fprintf(screen,"dz_Wji =  %.12f\n", dz_Wji);


        jmass = mass[jtype];

        ijmass = imass*jmass;

        rho_j2 = rho_SPH[j]*rho_SPH[j];

        // compute force terms from pressure tensor here
        
        f[i][0] += -ijmass*((Pixx*dx_Wij + Pixy*dy_Wij + Pixz*dz_Wij)/(rho_i2*omega_i) + (Pjxx*dx_Wji + Pjxy*dy_Wji + Pjxz*dz_Wji)/(rho_j2*omega_j));
        // fprintf(screen,"fx =  %.12f\n", -ijmass*((Pixx*dx_Wij + Pixy*dy_Wij + Pixz*dz_Wij)/(rho_i2*omega_i) + (Pjxx*dx_Wji + Pjxy*dy_Wji + Pjxz*dz_Wji)/(rho_j2*omega_j)));
        f[i][1] += -ijmass*((Pixy*dx_Wij + Piyy*dy_Wij + Piyz*dz_Wij)/(rho_i2*omega_i) + (Pjxy*dx_Wji + Pjyy*dy_Wji + Pjyz*dz_Wji)/(rho_j2*omega_j));
        // fprintf(screen,"fy =  %.12f\n", -ijmass*((Pixy*dx_Wij + Piyy*dy_Wij + Piyz*dz_Wij)/(rho_i2*omega_i) + (Pjxy*dx_Wji + Pjyy*dy_Wji + Pjyz*dz_Wji)/(rho_j2*omega_j)));
        f[i][2] += -ijmass*((Pixz*dx_Wij + Piyz*dy_Wij + Pizz*dz_Wij)/(rho_i2*omega_i) + (Pjxz*dx_Wji + Pjyz*dy_Wji + Pjzz*dz_Wji)/(rho_j2*omega_j));
        // fprintf(screen,"fz =  %.12f\n", -ijmass*((Pixz*dx_Wij + Piyz*dy_Wij + Pizz*dz_Wij)/(rho_i2*omega_i) + (Pjxz*dx_Wji + Pjyz*dy_Wji + Pjzz*dz_Wji)/(rho_j2*omega_j)));
        
        // fprintf(screen,"ijmass =  %.9f\n", ijmass);
        // fprintf(screen,"rho_j2 =  %.9f\n", rho_j2);
        // fprintf(screen,"omega_j =  %.9f\n", omega_j);
        
        // fprintf(screen,"dx_rho[j] =  %.9f\n", dx_rho[j]);
        // fprintf(screen,"dy_rho[j] =  %.9f\n", dy_rho[j]);
        // fprintf(screen,"dz_rho[j] =  %.9f\n", dz_rho[j]);
        // fprintf(screen,"dxx_rho[j] =  %.9f\n", dxx_rho[j]);
        // fprintf(screen,"dxy_rho[j] =  %.9f\n", dxy_rho[j]);
        // fprintf(screen,"dxz_rho[j] =  %.9f\n", dxz_rho[j]);
        // fprintf(screen,"dyz_rho[j] =  %.9f\n", dyz_rho[j]);
        // fprintf(screen,"dzz_rho[j] =  %.9f\n", dzz_rho[j]);

        // fprintf(screen,"Pjxx =  %.9f\n", Pjxx);
        // fprintf(screen,"Pjxy =  %.9f\n", Pjxy);
        // fprintf(screen,"Pjxz =  %.9f\n", Pjxz);
        // fprintf(screen,"Pjyy =  %.9f\n", Pjyy);
        // fprintf(screen,"Pjyz =  %.9f\n", Pjyz);
        // fprintf(screen,"Pjzz =  %.9f\n", Pjzz);

        // fprintf(screen,"dx_Wij =  %.9f\n", dx_Wij);
        // fprintf(screen,"dy_Wij =  %.9f\n", dy_Wij);
        // fprintf(screen,"dz_Wij =  %.9f\n", dz_Wij);

        // fprintf(screen,"dx_Wji =  %.9f\n", dx_Wji);
        // fprintf(screen,"dy_Wji =  %.9f\n", dy_Wji);
        // fprintf(screen,"dz_Wji =  %.9f\n", dz_Wji);
        
        if (newton_pair || j < nlocal) {
          // fprintf(screen,"break_statement H\n");
          f[j][0] += ijmass*((Pixx*dx_Wij + Pixy*dy_Wij + Pixz*dz_Wij)/(rho_i2*omega_i) + (Pjxx*dx_Wji + Pjxy*dy_Wji + Pjxz*dz_Wji)/(rho_j2*omega_j));
          f[j][1] += ijmass*((Pixy*dx_Wij + Piyy*dy_Wij + Piyz*dz_Wij)/(rho_i2*omega_i) + (Pjxy*dx_Wji + Pjyy*dy_Wji + Pjyz*dz_Wji)/(rho_j2*omega_j));
          f[j][2] += ijmass*((Pixz*dx_Wij + Piyz*dy_Wij + Pizz*dz_Wij)/(rho_i2*omega_i) + (Pjxz*dx_Wji + Pjyz*dy_Wji + Pjzz*dz_Wji)/(rho_j2*omega_j));
        
        }
      }
    }
  }
  if (vflag_fdotr) virial_fdotr_compute();
  // fprintf(screen,"break_statement I\n");
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairBohmSPHDynamicMocz::allocate()
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

void PairBohmSPHDynamicMocz::settings(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR,"Illegal pair_style command. Require 2 input arguments.");

  cut_global = force->numeric(FLERR,arg[0]);
  gamma_factor = force->numeric(FLERR,arg[1]);

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

void PairBohmSPHDynamicMocz::coeff(int narg, char **arg)
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

void PairBohmSPHDynamicMocz::init_style()
{
  neighbor->request(this,instance_me);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairBohmSPHDynamicMocz::init_one(int i, int j)
{
  if (setflag[i][j] == 0)
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);

  return cut[i][j];
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairBohmSPHDynamicMocz::write_restart(FILE *fp)
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

void PairBohmSPHDynamicMocz::read_restart(FILE *fp)
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

void PairBohmSPHDynamicMocz::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairBohmSPHDynamicMocz::read_restart_settings(FILE *fp)
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

int PairBohmSPHDynamicMocz::pack_forward_comm(int n, int *list, double *buf,
                               int /*pbc_flag*/, int * /*pbc*/)
{
  int i,j,m;

  m = 0;
  if (commflag == 0){
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = dx_rho[j];
      buf[m++] = dy_rho[j];
      buf[m++] = dz_rho[j];
      buf[m++] = dxx_rho[j];
      buf[m++] = dxy_rho[j];
      buf[m++] = dxz_rho[j];
      buf[m++] = dyy_rho[j];
      buf[m++] = dyz_rho[j];
      buf[m++] = dzz_rho[j];
    }
    return m;
  }
}

/* ---------------------------------------------------------------------- */

void PairBohmSPHDynamicMocz::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  if (commflag == 0){
    for (i = first; i < last; i++){
      dx_rho[i] = buf[m++];
      dy_rho[i] = buf[m++];
      dz_rho[i] = buf[m++];
      dxx_rho[i] = buf[m++];
      dxy_rho[i] = buf[m++];
      dxz_rho[i] = buf[m++];
      dyy_rho[i] = buf[m++];
      dyz_rho[i] = buf[m++];
      dzz_rho[i] = buf[m++];
    }
  }
}

/* ---------------------------------------------------------------------- */

int PairBohmSPHDynamicMocz::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  if (commflag == 0){
    for (i = first; i < last; i++){
      buf[m++] = dx_rho[i];
      buf[m++] = dy_rho[i];
      buf[m++] = dz_rho[i];
      buf[m++] = dxx_rho[i];
      buf[m++] = dxy_rho[i];
      buf[m++] = dxz_rho[i];
      buf[m++] = dyy_rho[i];
      buf[m++] = dyz_rho[i];
      buf[m++] = dzz_rho[i];
    }
    return m;
  }
}

/* ---------------------------------------------------------------------- */

void PairBohmSPHDynamicMocz::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  if (commflag == 0){
    for (i = 0; i < n; i++) {
      j = list[i];
      dx_rho[j] += buf[m++];
      dy_rho[j] += buf[m++];
      dz_rho[j] += buf[m++];
      dxx_rho[j] += buf[m++];
      dxy_rho[j] += buf[m++];
      dxz_rho[j] += buf[m++];
      dyy_rho[j] += buf[m++];
      dyz_rho[j] += buf[m++];
      dzz_rho[j] += buf[m++];
    }
  }
}

/* ---------------------------------------------------------------------- */

void *PairBohmSPHDynamicMocz::extract(const char *str, int &dim)
{
  dim = 2;
  return NULL;
}
