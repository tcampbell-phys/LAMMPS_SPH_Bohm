/* -*- c++ -*- ----------------------------------------------------------
Bespoke pair_style to compute the Bohm force using an SPH-style pressure
tensor.

pair_bohm_SPH_basic:
- Electrons as basic fluid element.
- Static, global, gaussian widths.

Thomas Campbell
------------------------------------------------------------------------- */
#include "pair_bohm_SPH_basic.h"
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

PairBohmSPHBasic::PairBohmSPHBasic(LAMMPS *lmp) : Pair(lmp) {
  nmax = 0;

  manybody_flag = 1;

  rho = NULL;
  dx_rho = NULL;
  dy_rho = NULL;
  dz_rho = NULL;
  dxx_rho = NULL;
  dxy_rho = NULL;
  dxz_rho = NULL;
  dyy_rho = NULL;
  dyz_rho = NULL;
  dzz_rho = NULL;

  comm_forward = 10;
  comm_reverse = 10;
}

/* ---------------------------------------------------------------------- */

PairBohmSPHBasic::~PairBohmSPHBasic()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
    memory->destroy(rho);
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

void PairBohmSPHBasic::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz;
  double xtmp2,ytmp2,ztmp2,delx2,dely2,delz2;
  double delx_2,dely_2,delz_2,rsq;
  double gauss_pre;
  double h2,hm2,hm4;
  double imass,jmass,ijmass;
  double m_gauss_ij,m_gauss_ji;
  double omega_ij;
  double Pixx,Pixy,Pixz,Piyy,Piyz,Pizz;
  double Pjxx,Pjxy,Pjxz,Pjyy,Pjyz,Pjzz;
  double dx_Wij,dy_Wij,dz_Wij;
  double rho_i2,rho_j2;
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
    memory->destroy(rho);
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

    memory->create(rho,nmax,"pair:rho");
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
      rho[i] = 0.0;
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
      rho[i] = 0.0;
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

  // 3D Gaussian prefactor
  gauss_pre = 1/(pow(2*M_PI,1.5)*pow(width_global,3));

  // width terms
  h2 = width_global*width_global;
  hm2 = 1/h2;
  hm4 = hm2*hm2;

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

    // self interaction terms
    rho[i] += imass*gauss_pre;
    dxx_rho[i] += -hm2*imass*gauss_pre;
    dyy_rho[i] += -hm2*imass*gauss_pre;
    dzz_rho[i] += -hm2*imass*gauss_pre;

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

        m_gauss_ij = jmass*gauss_pre*exp(-(rsq)/(2*h2));
      
        rho[i] += m_gauss_ij;
        dx_rho[i] += ((-delx)*hm2)*m_gauss_ij;
        dy_rho[i] += ((-dely)*hm2)*m_gauss_ij;
        dz_rho[i] += ((-delz)*hm2)*m_gauss_ij;
        dxx_rho[i] += hm2*(delx_2*hm2 - 1)*m_gauss_ij;
        dxy_rho[i] += hm4*(delx*dely)*m_gauss_ij;
        dxz_rho[i] += hm4*(delx*delz)*m_gauss_ij;
        dyy_rho[i] += hm2*(dely_2*hm2 - 1)*m_gauss_ij;
        dyz_rho[i] += hm4*(dely*delz)*m_gauss_ij;
        dzz_rho[i] += hm2*(delz_2*hm2 - 1)*m_gauss_ij;

        if (newton_pair || j < nlocal) {
        
          m_gauss_ji = imass*gauss_pre*exp(-(rsq)/(2*h2));

          rho[j] += m_gauss_ji;
          dx_rho[j] += ((delx)*hm2)*m_gauss_ji;
          dy_rho[j] += ((dely)*hm2)*m_gauss_ji;
          dz_rho[j] += ((delz)*hm2)*m_gauss_ji;
          dxx_rho[j] += hm2*(delx_2*hm2 - 1)*m_gauss_ji;
          dxy_rho[j] += hm4*(delx*dely)*m_gauss_ji;
          dxz_rho[j] += hm4*(delx*delz)*m_gauss_ji;
          dyy_rho[j] += hm2*(dely_2*hm2 - 1)*m_gauss_ji;
          dyz_rho[j] += hm4*(dely*delz)*m_gauss_ji;
          dzz_rho[j] += hm2*(delz_2*hm2 - 1)*m_gauss_ji;
        }
      }
    }
  }
  
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

    rho_i2 = rho[i]*rho[i];

    Pixx = gamma_factor*f_prefactor*((dx_rho[i]*dx_rho[i])/rho[i] - dxx_rho[i]);
    Pixy = gamma_factor*f_prefactor*((dx_rho[i]*dy_rho[i])/rho[i] - dxy_rho[i]);
    Pixz = gamma_factor*f_prefactor*((dx_rho[i]*dz_rho[i])/rho[i] - dxz_rho[i]);
    Piyy = gamma_factor*f_prefactor*((dy_rho[i]*dy_rho[i])/rho[i] - dyy_rho[i]);
    Piyz = gamma_factor*f_prefactor*((dy_rho[i]*dz_rho[i])/rho[i] - dyz_rho[i]);
    Pizz = gamma_factor*f_prefactor*((dz_rho[i]*dz_rho[i])/rho[i] - dzz_rho[i]);

    // Bohm potential calculation
    bohm_pot = -gamma_factor*f_prefactor*((dxx_rho[i] + dyy_rho[i]+ dzz_rho[i])/rho[i] - (dx_rho[i]*dx_rho[i] + dy_rho[i]*dy_rho[i] + dz_rho[i]*dz_rho[i])/(2*rho_i2));
    if (eflag_global) eng_vdwl += bohm_pot;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];

      rsq = delx*delx + dely*dely + delz*delz;

      jtype = type[j];
 
      if (rsq < cutsquared) {

        omega_ij = exp(-(rsq)/(2*h2));

        Pjxx = gamma_factor*f_prefactor*((dx_rho[j]*dx_rho[j])/rho[j] - dxx_rho[j]);
        Pjxy = gamma_factor*f_prefactor*((dx_rho[j]*dy_rho[j])/rho[j] - dxy_rho[j]);
        Pjxz = gamma_factor*f_prefactor*((dx_rho[j]*dz_rho[j])/rho[j] - dxz_rho[j]);
        Pjyy = gamma_factor*f_prefactor*((dy_rho[j]*dy_rho[j])/rho[j] - dyy_rho[j]);
        Pjyz = gamma_factor*f_prefactor*((dy_rho[j]*dz_rho[j])/rho[j] - dyz_rho[j]);
        Pjzz = gamma_factor*f_prefactor*((dz_rho[j]*dz_rho[j])/rho[j] - dzz_rho[j]);
      
        dx_Wij = gauss_pre*(-delx/h2)*omega_ij;
        dy_Wij = gauss_pre*(-dely/h2)*omega_ij;
        dz_Wij = gauss_pre*(-delz/h2)*omega_ij;

        jmass = mass[jtype];

        ijmass = imass*jmass;

        rho_j2 = rho[j]*rho[j];

        // compute force terms from pressure tensor here
        
        f[i][0] += -ijmass*((Pixx*dx_Wij + Pixy*dy_Wij + Pixz*dz_Wij)/rho_i2 + (Pjxx*dx_Wij + Pjxy*dy_Wij + Pjxz*dz_Wij)/rho_j2);
        f[i][1] += -ijmass*((Pixy*dx_Wij + Piyy*dy_Wij + Piyz*dz_Wij)/rho_i2 + (Pjxy*dx_Wij + Pjyy*dy_Wij + Pjyz*dz_Wij)/rho_j2);
        f[i][2] += -ijmass*((Pixz*dx_Wij + Piyz*dy_Wij + Pizz*dz_Wij)/rho_i2 + (Pjxz*dx_Wij + Pjyz*dy_Wij + Pjzz*dz_Wij)/rho_j2);
        
        if (newton_pair || j < nlocal) {
          f[j][0] += ijmass*((Pixx*dx_Wij + Pixy*dy_Wij + Pixz*dz_Wij)/rho_i2 + (Pjxx*dx_Wij + Pjxy*dy_Wij + Pjxz*dz_Wij)/rho_j2);
          f[j][1] += ijmass*((Pixy*dx_Wij + Piyy*dy_Wij + Piyz*dz_Wij)/rho_i2 + (Pjxy*dx_Wij + Pjyy*dy_Wij + Pjyz*dz_Wij)/rho_j2);
          f[j][2] += ijmass*((Pixz*dx_Wij + Piyz*dy_Wij + Pizz*dz_Wij)/rho_i2 + (Pjxz*dx_Wij + Pjyz*dy_Wij + Pjzz*dz_Wij)/rho_j2);
          
        }
      }
    }
  }
  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairBohmSPHBasic::allocate()
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

void PairBohmSPHBasic::settings(int narg, char **arg)
{
  if (narg != 3) error->all(FLERR,"Illegal pair_style command. Require 3 input arguments.");

  cut_global = force->numeric(FLERR,arg[0]);
  width_global = force->numeric(FLERR,arg[1]);
  gamma_factor = force->numeric(FLERR,arg[2]);

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

void PairBohmSPHBasic::coeff(int narg, char **arg)
{
  if (narg < 2 || narg > 4)
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

// Specify a different sigma in pair_coeff step:

  double sigma_one = width_global;
  if (narg == 4) sigma_one = force->numeric(FLERR,arg[3]);

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}


/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairBohmSPHBasic::init_style()
{
  neighbor->request(this,instance_me);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairBohmSPHBasic::init_one(int i, int j)
{
  if (setflag[i][j] == 0)
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);

  return cut[i][j];
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairBohmSPHBasic::write_restart(FILE *fp)
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

void PairBohmSPHBasic::read_restart(FILE *fp)
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

void PairBohmSPHBasic::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairBohmSPHBasic::read_restart_settings(FILE *fp)
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

int PairBohmSPHBasic::pack_forward_comm(int n, int *list, double *buf,
                               int /*pbc_flag*/, int * /*pbc*/)
{
  int i,j,m;

  m = 0;
  if (commflag == 0){
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = rho[j];
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

void PairBohmSPHBasic::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  if (commflag == 0){
    for (i = first; i < last; i++){
      rho[i] = buf[m++];
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

int PairBohmSPHBasic::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  if (commflag == 0){
    for (i = first; i < last; i++){
      buf[m++] = rho[i];
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

void PairBohmSPHBasic::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  if (commflag == 0){
    for (i = 0; i < n; i++) {
      j = list[i];
      rho[j] += buf[m++];
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

void *PairBohmSPHBasic::extract(const char *str, int &dim)
{
  dim = 2;
  return NULL;
}
