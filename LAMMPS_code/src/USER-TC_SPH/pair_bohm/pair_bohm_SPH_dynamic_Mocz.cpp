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

  dxx_rho = NULL;
  dxy_rho = NULL;
  dxz_rho = NULL;
  dyy_rho = NULL;
  dyz_rho = NULL;
  dzz_rho = NULL;

  comm_forward = 6;
  comm_reverse = 6;
}

/* ---------------------------------------------------------------------- */

PairBohmSPHDynamicMocz::~PairBohmSPHDynamicMocz()
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
  }
}

/* ---------------------------------------------------------------------- */

void PairBohmSPHDynamicMocz::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz;
  double vxtmp,vytmp,vztmp,delvx,delvy,delvz;
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
  double u_prefact_i,u_prefact_j;

  ev_init(eflag,vflag);
  
  hplanck  = force->hplanck;
  hbar = hplanck/(2*M_PI);

  // Bohm pressure prefactor
  f_prefactor = (hbar*hbar)/(4*e_mass);

  int *ilist,*jlist,*numneigh,**firstneigh;

  if (atom->nmax > nmax) {
    // delete and create new memory arrays for any per-particle variables that need communicating.
    memory->destroy(dxx_rho);
    memory->destroy(dxy_rho);
    memory->destroy(dxz_rho);
    memory->destroy(dyy_rho);
    memory->destroy(dyz_rho);
    memory->destroy(dzz_rho);
    
    nmax = atom->nmax;

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
  double *dx_rho_SPH = atom->dx_rho_SPH;
  double *dy_rho_SPH = atom->dy_rho_SPH;
  double *dz_rho_SPH = atom->dz_rho_SPH;
  double *u_SPH = atom->u_SPH;

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

  cutsquared = cut_global*cut_global;

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
        
        dxx_rho[i] += hm2_i*(delx_2*hm2_i - 1.)*m_gauss_ij*dens_factor_ij;
        dxy_rho[i] += hm4_i*(delx*dely)*m_gauss_ij*dens_factor_ij;
        dxz_rho[i] += hm4_i*(delx*delz)*m_gauss_ij*dens_factor_ij;
        dyy_rho[i] += hm2_i*(dely_2*hm2_i - 1.)*m_gauss_ij*dens_factor_ij;
        dyz_rho[i] += hm4_i*(dely*delz)*m_gauss_ij*dens_factor_ij;
        dzz_rho[i] += hm2_i*(delz_2*hm2_i - 1.)*m_gauss_ij*dens_factor_ij;

        if (newton_pair || j < nlocal) {

          h_j = width_SPH[j];
          h2_j = h_j*h_j;
          hm2_j = 1./h2_j;
          hm4_j = hm2_j*hm2_j;

          dens_factor_ji = 1. - rho_j/rho_i;

          gauss_pre_j = pi_fact*(1./(h_j*h_j*h_j));
        
          m_gauss_ji = imass*gauss_pre_j*exp(-(rsq)*hm2_j/2.);

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
  comm_forward = 6;
  comm_reverse = 6;
  commflag = 0;

  if (newton_pair) comm->reverse_comm_pair(this);

  comm->forward_comm_pair(this);

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

    // // fprintf(screen,"\nParticle at...\n");
    // // fprintf(screen,"x = %16.16f\n",xtmp);
    // // fprintf(screen,"y = %16.16f\n",ytmp);
    // // fprintf(screen,"z = %16.16f\n",ztmp);

    // // fprintf(screen,"dx_rho_SPH = %16.16f\n",dx_rho_SPH[i]);
    // // fprintf(screen,"dy_rho_SPH = %16.16f\n",dy_rho_SPH[i]);
    // // fprintf(screen,"dz_rho_SPH = %16.16f\n",dz_rho_SPH[i]);

    Pixx = gamma_factor*f_prefactor*((dx_rho_SPH[i]*dx_rho_SPH[i])/rho_SPH[i] - dxx_rho[i]);
    Pixy = gamma_factor*f_prefactor*((dx_rho_SPH[i]*dy_rho_SPH[i])/rho_SPH[i] - dxy_rho[i]);
    Pixz = gamma_factor*f_prefactor*((dx_rho_SPH[i]*dz_rho_SPH[i])/rho_SPH[i] - dxz_rho[i]);
    Piyy = gamma_factor*f_prefactor*((dy_rho_SPH[i]*dy_rho_SPH[i])/rho_SPH[i] - dyy_rho[i]);
    Piyz = gamma_factor*f_prefactor*((dy_rho_SPH[i]*dz_rho_SPH[i])/rho_SPH[i] - dyz_rho[i]);
    Pizz = gamma_factor*f_prefactor*((dz_rho_SPH[i]*dz_rho_SPH[i])/rho_SPH[i] - dzz_rho[i]);

    u_prefact_i = (dt/(rho_i2*omega_i));

    // fprintf(screen,"\nBohm h_i = %16.16f",h_i);

    // fprintf(screen,"\nBohm omega_i = %16.16f",omega_i);

    // fprintf(screen,"\nBohm rho_i2 = %16.16f",rho_i2);

    // fprintf(screen,"\nBohm Pixx = %16.16f",Pixx);

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

        Pjxx = gamma_factor*f_prefactor*((dx_rho_SPH[j]*dx_rho_SPH[j])/rho_SPH[j] - dxx_rho[j]);
        Pjxy = gamma_factor*f_prefactor*((dx_rho_SPH[j]*dy_rho_SPH[j])/rho_SPH[j] - dxy_rho[j]);
        Pjxz = gamma_factor*f_prefactor*((dx_rho_SPH[j]*dz_rho_SPH[j])/rho_SPH[j] - dxz_rho[j]);
        Pjyy = gamma_factor*f_prefactor*((dy_rho_SPH[j]*dy_rho_SPH[j])/rho_SPH[j] - dyy_rho[j]);
        Pjyz = gamma_factor*f_prefactor*((dy_rho_SPH[j]*dz_rho_SPH[j])/rho_SPH[j] - dyz_rho[j]);
        Pjzz = gamma_factor*f_prefactor*((dz_rho_SPH[j]*dz_rho_SPH[j])/rho_SPH[j] - dzz_rho[j]);
      
        dx_Wij = gauss_pre_i*(-delx*hm2_i)*exp_ij;
        dy_Wij = gauss_pre_i*(-dely*hm2_i)*exp_ij;
        dz_Wij = gauss_pre_i*(-delz*hm2_i)*exp_ij;

        dx_Wji = gauss_pre_j*(-delx*hm2_j)*exp_ji;
        dy_Wji = gauss_pre_j*(-dely*hm2_j)*exp_ji;
        dz_Wji = gauss_pre_j*(-delz*hm2_j)*exp_ji;

        rho_j2 = rho_SPH[j]*rho_SPH[j];

        u_prefact_j = (dt/(rho_j2*omega_j));

        u_SPH[i] += u_prefact_i*ijmass*((Pixx*delvx + Pixy*delvy + Pixz*delvz)*dx_Wij + (Pixy*delvx + Piyy*delvy + Piyz*delvz)*dy_Wij + (Pixz*delvx + Piyz*delvy + Pizz*delvz)*dz_Wij);
        
        // compute force terms from pressure tensor here

        fx = -ijmass*((Pixx*dx_Wij + Pixy*dy_Wij + Pixz*dz_Wij)/(rho_i2*omega_i) + (Pjxx*dx_Wji + Pjxy*dy_Wji + Pjxz*dz_Wji)/(rho_j2*omega_j));
        fy = -ijmass*((Pixy*dx_Wij + Piyy*dy_Wij + Piyz*dz_Wij)/(rho_i2*omega_i) + (Pjxy*dx_Wji + Pjyy*dy_Wji + Pjyz*dz_Wji)/(rho_j2*omega_j));
        fz = -ijmass*((Pixz*dx_Wij + Piyz*dy_Wij + Pizz*dz_Wij)/(rho_i2*omega_i) + (Pjxz*dx_Wji + Pjyz*dy_Wji + Pjzz*dz_Wji)/(rho_j2*omega_j));
        
        f[i][0] += fx;
        f[i][1] += fy;
        f[i][2] += fz;

        // fprintf(screen,"\nBohm fx = %16.16f",fx);

        if (newton_pair || j < nlocal) {
          f[j][0] -= fx;
          f[j][1] -= fy;
          f[j][2] -= fz;

          u_SPH[j] -= u_prefact_j*ijmass*((-Pjxx*delvx - Pjxy*delvy - Pjxz*delvz)*dx_Wji + (-Pjxy*delvx - Pjyy*delvy - Pjyz*delvz)*dy_Wji + (-Pjxz*delvx - Pjyz*delvy - Pjzz*delvz)*dz_Wji);
        }

        eflag_either = 0;
        if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,0.0,0.0,
                        fx,fy,fz,delx,dely,delz);
      }
    }
  }
  comm_forward = 1;
  comm_reverse = 1;
  commflag = 1;
  if (newton_pair) comm->reverse_comm_pair(this);
  comm->forward_comm_pair(this);

  eflag_either = 1;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (eflag_global) eng_vdwl +=  u_SPH[i];
    if (eflag_atom) eatom[i] += u_SPH[i];
  }

  if (vflag_fdotr) virial_fdotr_compute();
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

  if (comm->ghost_velocity != 1) error->all(FLERR,"Illegal pair_style command. Require ghost atom velocity.");

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
      buf[m++] = dxx_rho[j];
      buf[m++] = dxy_rho[j];
      buf[m++] = dxz_rho[j];
      buf[m++] = dyy_rho[j];
      buf[m++] = dyz_rho[j];
      buf[m++] = dzz_rho[j];
    }
    return m;
  }
  if (commflag == 1){
    double *u_SPH = atom->u_SPH;
    // // fprintf(screen,"In u_SPH pack forward comm loop...");
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = u_SPH[j];
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
      dxx_rho[i] = buf[m++];
      dxy_rho[i] = buf[m++];
      dxz_rho[i] = buf[m++];
      dyy_rho[i] = buf[m++];
      dyz_rho[i] = buf[m++];
      dzz_rho[i] = buf[m++];
    }
  }
  if (commflag == 1){
    double *u_SPH = atom->u_SPH;
    // // fprintf(screen,"In u_SPH unpack forward comm loop...");
    for (i = first; i < last; i++){
      u_SPH[i] = buf[m++];
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
      buf[m++] = dxx_rho[i];
      buf[m++] = dxy_rho[i];
      buf[m++] = dxz_rho[i];
      buf[m++] = dyy_rho[i];
      buf[m++] = dyz_rho[i];
      buf[m++] = dzz_rho[i];
    }
    return m;
  }
  if (commflag == 1){
    double *u_SPH = atom->u_SPH;
    // // fprintf(screen,"In u_SPH pack reverse comm loop...");
    for (i = first; i < last; i++){
      buf[m++] = u_SPH[i];
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
      dxx_rho[j] += buf[m++];
      dxy_rho[j] += buf[m++];
      dxz_rho[j] += buf[m++];
      dyy_rho[j] += buf[m++];
      dyz_rho[j] += buf[m++];
      dzz_rho[j] += buf[m++];
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
}

/* ---------------------------------------------------------------------- */

void *PairBohmSPHDynamicMocz::extract(const char *str, int &dim)
{
  dim = 2;
  return NULL;
}
