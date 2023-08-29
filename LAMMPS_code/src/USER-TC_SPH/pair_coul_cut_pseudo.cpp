/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "pair_coul_cut_pseudo.h"
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

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairCoulCutPseudo::PairCoulCutPseudo(LAMMPS *lmp) : Pair(lmp) {
  centroidstressflag = 1;
}

/* ---------------------------------------------------------------------- */

PairCoulCutPseudo::~PairCoulCutPseudo()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(scale);
  }
}

/* ---------------------------------------------------------------------- */

void PairCoulCutPseudo::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype,k;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,ecoul,fpair;
  double rsq,r2inv,rinv,forcecoul,factor_coul;
  int *ilist,*jlist,*numneigh,**firstneigh;

  ecoul = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;

  double force_fact;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  //// fprintf(screen,"qqrd2e = %f \n",qqrd2e);

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

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_coul = special_coul[sbmask(j)];
      //// fprintf(screen,"factor_coul = %f \n",factor_coul);
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        r2inv = 1.0/rsq;
        rinv = sqrt(r2inv);

        if (jtype != ion_species){

          // electron neighbour
          if (itype == ion_species){
            // fprintf(screen,"\nrsq = %16.16f",rsq);
            // fprintf(screen,"\ndx = %16.16f",delx);
            // ion target
            ecoul = 0.;
            force_fact = 0.;
            for (k = 0; k < N_coeff; k++){
              ecoul += c_coeff[k] * exp(-a_coeff[k]*rsq);
              force_fact += 2 * c_coeff[k] * a_coeff[k] * exp(-a_coeff[k]*rsq);
            }
            ecoul *= factor_coul * qqrd2e * scale[itype][jtype];
            force_fact *= factor_coul * qqrd2e * scale[itype][jtype];
            // fprintf(screen,"\necoul = %16.16f",ecoul);
            // fprintf(screen,"\nfx = %16.16f",force_fact*delx);
            f[i][0] += delx*force_fact;
            f[i][1] += dely*force_fact;
            f[i][2] += delz*force_fact;

            if (newton_pair || j < nlocal) {
              f[j][0] -= delx*force_fact;
              f[j][1] -= dely*force_fact;
              f[j][2] -= delz*force_fact;

            }
            if (evflag) ev_tally(i,j,nlocal,newton_pair,
                                 0.0,ecoul,force_fact,delx,dely,delz);
          }
        }
        if (jtype == ion_species){

          // ion neighbour

          if (itype != ion_species){
            // fprintf(screen,"\nrsq = %16.16f",rsq);
            // fprintf(screen,"\ndx = %16.16f",delx);

            // electron target
            ecoul = 0.;
            force_fact = 0.;

            for (k = 0; k < N_coeff; k++){
              ecoul += c_coeff[k] * exp(-a_coeff[k]*rsq);
              force_fact += 2 * c_coeff[k] * a_coeff[k] * exp(-a_coeff[k]*rsq);
            }
            
            ecoul *= factor_coul * qqrd2e * scale[itype][jtype];
            force_fact *= factor_coul * qqrd2e * scale[itype][jtype];
            // fprintf(screen,"\necoul = %16.16f",ecoul);
            // fprintf(screen,"\nfx = %16.16f",force_fact*delx);
            f[i][0] += delx*force_fact;
            f[i][1] += dely*force_fact;
            f[i][2] += delz*force_fact;

            if (newton_pair || j < nlocal) {
              f[j][0] -= delx*force_fact;
              f[j][1] -= dely*force_fact;
              f[j][2] -= delz*force_fact;

            }
            if (evflag) ev_tally(i,j,nlocal,newton_pair,
                                 0.0,ecoul,force_fact,delx,dely,delz);
          }
        }
        if (jtype == ion_species){

          // ion neighbour

          if (itype == ion_species){

            // ion target
            forcecoul = qqrd2e * scale[itype][jtype] * qtmp*q[j]*rinv;
            //// fprintf(screen,"scale[itype][jtype] = %f \n",scale[itype][jtype]);
            fpair = factor_coul*forcecoul * r2inv;

            f[i][0] += delx*fpair;
            f[i][1] += dely*fpair;
            f[i][2] += delz*fpair;
            if (newton_pair || j < nlocal) {
              f[j][0] -= delx*fpair;
              f[j][1] -= dely*fpair;
              f[j][2] -= delz*fpair;
            }

            if (eflag)
              ecoul = factor_coul * qqrd2e * scale[itype][jtype] * qtmp*q[j]*rinv;

            if (evflag) ev_tally(i,j,nlocal,newton_pair,
                                0.0,ecoul,fpair,delx,dely,delz);
          }
        }
        if (jtype != ion_species){

          // electron neighbour

          if (itype != ion_species){

            // electron target
            forcecoul = qqrd2e * scale[itype][jtype] * qtmp*q[j]*rinv;
            //// fprintf(screen,"scale[itype][jtype] = %f \n",scale[itype][jtype]);
            fpair = factor_coul*forcecoul * r2inv;

            f[i][0] += delx*fpair;
            f[i][1] += dely*fpair;
            f[i][2] += delz*fpair;
            if (newton_pair || j < nlocal) {
              f[j][0] -= delx*fpair;
              f[j][1] -= dely*fpair;
              f[j][2] -= delz*fpair;
            }

            if (eflag)
              ecoul = factor_coul * qqrd2e * scale[itype][jtype] * qtmp*q[j]*rinv;

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

void PairCoulCutPseudo::allocate()
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

void PairCoulCutPseudo::settings(int narg, char **arg)
{
  if (narg != 4) error->all(FLERR,"Illegal pair_style command");

  cut_global = force->numeric(FLERR,arg[0]);
  ion_species = force->numeric(FLERR,arg[1]);
  pseudo_key = force->numeric(FLERR,arg[2]);
  ion_charge = force->numeric(FLERR,arg[3]);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
  if (pseudo_key==al_lda_D_key){
    c_coeff = al_lda_D_c;
    a_coeff = al_lda_D_a;
    N_coeff = al_lda_D_Ncoeff;
    if (ion_charge != al_lda_D_ion_charge){
      // fprintf(screen,"\n ion_charge = %d \n",ion_charge);
      // fprintf(screen,"\n al_lda_D_ion_charge = %d \n",al_lda_D_ion_charge);
      error->all(FLERR,"PSEUDO_ERR requested pseudopotential parameters do not match input: ion_charge ");
    }
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairCoulCutPseudo::coeff(int narg, char **arg)
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

void PairCoulCutPseudo::init_style()
{
  if (!atom->q_flag)
    error->all(FLERR,"Pair style coul/cut requires atom attribute q");

  neighbor->request(this,instance_me);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairCoulCutPseudo::init_one(int i, int j)
{
  if (setflag[i][j] == 0)
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);

  scale[j][i] = scale[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairCoulCutPseudo::write_restart(FILE *fp)
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

void PairCoulCutPseudo::read_restart(FILE *fp)
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

void PairCoulCutPseudo::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairCoulCutPseudo::read_restart_settings(FILE *fp)
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

double PairCoulCutPseudo::single(int i, int j, int /*itype*/, int /*jtype*/,
                           double rsq, double factor_coul, double /*factor_lj*/,
                           double &fforce)
{
  double r2inv,rinv,forcecoul,phicoul;

  r2inv = 1.0/rsq;
  rinv = sqrt(r2inv);
  forcecoul = force->qqrd2e * atom->q[i]*atom->q[j]*rinv;
  fforce = factor_coul*forcecoul * r2inv;

  phicoul = force->qqrd2e * atom->q[i]*atom->q[j]*rinv;
  return factor_coul*phicoul;
}

/* ---------------------------------------------------------------------- */

void *PairCoulCutPseudo::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str,"cut_coul") == 0) return (void *) &cut;
  if (strcmp(str,"scale") == 0) return (void *) scale;
  return NULL;
}
