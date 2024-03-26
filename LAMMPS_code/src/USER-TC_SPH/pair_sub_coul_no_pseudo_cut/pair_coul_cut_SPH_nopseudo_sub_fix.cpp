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

#include "pair_coul_cut_SPH_nopseudo_sub_fix.h"
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
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairCoulCutSPHNoPseudoSubFix::PairCoulCutSPHNoPseudoSubFix(LAMMPS *lmp) : Pair(lmp) {
  nmax = 0;

  centroidstressflag = 1;
  manybody_flag = 1;

  comm_forward = 2;
  comm_reverse = 2;

}

/* ---------------------------------------------------------------------- */

PairCoulCutSPHNoPseudoSubFix::~PairCoulCutSPHNoPseudoSubFix()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(scale);
  }
}

/* ---------------------------------------------------------------------- */

void PairCoulCutSPHNoPseudoSubFix::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,ecoul,fpair;
  double rsq,rsqrt,r2inv,rinv,forcecoul,factor_coul;
  int *ilist,*jlist,*numneigh,**firstneigh;

  ecoul = 0.0;
  ev_init(eflag,vflag);

  int *tagid = atom->tag;
  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  double *special_coul = force->special_coul;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;
  double *width_SPH = atom->width_SPH;

  double rho_i, omega_i, fact_i;

  double eff_width;
  double force_fact;

  double h_i,h2_i,hm2_i,hm4_i;
  double h_j,h2_j,hm2_j,hm4_j;
  double imass,jmass,ijmass;
  double full_factor;
  double *mass = atom->mass;

  int ntimestep = update->ntimestep; // current timestep

  int lo_lim_lev;
  int hi_lim_lev;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

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

      // identify self-interaction particles
      lo_lim_lev = floor((tagid[i]-tag_ele_start)/N_epe)*N_epe + tag_ele_start;
      hi_lim_lev = lo_lim_lev + N_epe;

      h_i = width_SPH[i];
      h2_i = h_i*h_i;
      hm2_i = 1/h2_i;
   
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

            full_factor = 0.0;

            ecoul = 0.0;

            if ( tagid[j] >= hi_lim_lev || tagid[j] < lo_lim_lev ){

              eff_width = pow((h2_i + width_SPH[j]*width_SPH[j]),0.5);

              force_fact = factor_coul*qqrd2e*scale[itype][jtype]*(qtmp*q[j])*(erf(rsqrt/(sqrt2*eff_width))/(rsqrt*rsqrt*rsqrt) - (sqrt2/sqrt_pi)*(exp(-rsq/(2*eff_width*eff_width))/(eff_width*rsq)));

              f[i][0] += delx*force_fact;
              f[i][1] += dely*force_fact;
              f[i][2] += delz*force_fact;

              if (newton_pair || j < nlocal) {

                f[j][0] -= delx*force_fact;
                f[j][1] -= dely*force_fact;
                f[j][2] -= delz*force_fact;
              }

              full_factor = force_fact;

              ecoul = factor_coul * qqrd2e*scale[itype][jtype]*qtmp*q[j]*erf(rsqrt/(sqrt2*eff_width))/rsqrt;

            }

            // dynamic coulomb-SPH force expression is not pairwise symmetric
            // use of ev_tally not accurate for pressure evaluation - edit in future.
            if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,0.0,ecoul,
                        full_factor*delx,full_factor*dely,full_factor*delz,delx,dely,delz);
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

void PairCoulCutSPHNoPseudoSubFix::allocate()
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

void PairCoulCutSPHNoPseudoSubFix::settings(int narg, char **arg)
{
  if (narg != 7) error->all(FLERR,"Illegal pair_style command");

  cut_global = force->numeric(FLERR,arg[0]);
  ke_in = force->numeric(FLERR,arg[1]);
  ion_species = force->numeric(FLERR,arg[2]);
  N_epe = force->numeric(FLERR,arg[3]);
  tag_ele_start = force->numeric(FLERR,arg[4]);
  assign_ion = force->numeric(FLERR,arg[5]);
  targ_coord = force->numeric(FLERR,arg[6]);
  

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

void PairCoulCutSPHNoPseudoSubFix::coeff(int narg, char **arg)
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

void PairCoulCutSPHNoPseudoSubFix::init_style()
{
  if (!atom->q_flag)
    error->all(FLERR,"Pair style coul/cut requires atom attribute q");

  neighbor->request(this,instance_me);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairCoulCutSPHNoPseudoSubFix::init_one(int i, int j)
{
  if (setflag[i][j] == 0)
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);

  scale[j][i] = scale[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairCoulCutSPHNoPseudoSubFix::write_restart(FILE *fp)
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

void PairCoulCutSPHNoPseudoSubFix::read_restart(FILE *fp)
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

void PairCoulCutSPHNoPseudoSubFix::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairCoulCutSPHNoPseudoSubFix::read_restart_settings(FILE *fp)
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

void *PairCoulCutSPHNoPseudoSubFix::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str,"cut_coul") == 0) return (void *) &cut;
  if (strcmp(str,"scale") == 0) return (void *) scale;
  return NULL;
}

/* ---------------------------------------------------------------------- */


int PairCoulCutSPHNoPseudoSubFix::pack_forward_comm(int n, int *list, double *buf,
                               int /*pbc_flag*/, int * /*pbc*/)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = theta_coul[j];
    buf[m++] = theta_coul_ei[j];
  }
  return m;
  
}

/* ---------------------------------------------------------------------- */

void PairCoulCutSPHNoPseudoSubFix::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++){
    theta_coul[i] = buf[m++];
    theta_coul_ei[i] = buf[m++];
  }
  
}
/* ---------------------------------------------------------------------- */

int PairCoulCutSPHNoPseudoSubFix::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++){
    buf[m++] = theta_coul[i];
    buf[m++] = theta_coul_ei[i];
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void PairCoulCutSPHNoPseudoSubFix::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    theta_coul[j] += buf[m++];
    theta_coul_ei[j] += buf[m++];
  }
  
}
