/* -*- c++ -*- ----------------------------------------------------------------------
Thomas Campbell (Oxford)
------------------------------------------------------------------------- */

#include "pair_confine.h"
#include "domain.h"
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

PairConfine::PairConfine(LAMMPS *lmp) : Pair(lmp)
{
  boltz_val = force->boltz;
  hbar_val = (force->hplanck)/(2*M_PI);
  // fprintf(screen,"\nIn pair_confine with correct forcing term...");
}

/* ---------------------------------------------------------------------- */

PairConfine::~PairConfine()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
  }
}

/* ---------------------------------------------------------------------- */

void PairConfine::compute(int eflag, int vflag)
{
  int i,ii,inum,itype;
  double delx,dely,delz;
  int *ilist;

  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;

  double e_confine;

  int *tagid = atom->tag;
  // int ntimestep = update->ntimestep;

  double *x_COM = atom->x_COM;
  double *y_COM = atom->y_COM;
  double *z_COM = atom->z_COM;

  double x_COM_use,y_COM_use,z_COM_use;

  double add_factor_x,add_factor_y,add_factor_z;

  inum = list->inum;
  ilist = list->ilist;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];

    // Choose correct projection of Centre of Mass

    if (x[i][0]-x_COM[i] > half_box_len){
      // fprintf(screen,"\nBasic Ax");
      x_COM_use = x_COM[i] + box_len;
    } else if (x[i][0]-x_COM[i] < -half_box_len){
      // fprintf(screen,"\nBasic Bx");
      x_COM_use = x_COM[i] - box_len;
    } else{
      x_COM_use = x_COM[i];
    }
    if (x[i][1]-y_COM[i] > half_box_len){
      // fprintf(screen,"\nBasic Ay");
      y_COM_use = y_COM[i] + box_len;
    } else if (x[i][1]-y_COM[i] < -half_box_len){
      // fprintf(screen,"\nBasic By");
      y_COM_use = y_COM[i] - box_len;
    } else{
      y_COM_use = y_COM[i];
    }
    if (x[i][2]-z_COM[i] > half_box_len){
      // fprintf(screen,"\nBasic Az");
      z_COM_use = z_COM[i] + box_len;
    } else if (x[i][2]-z_COM[i] < -half_box_len){
      // fprintf(screen,"\nBasic Bz");
      z_COM_use = z_COM[i] - box_len;
    } else{
      z_COM_use = z_COM[i];
    }

    // if (ntimestep < 10){
    //   if (tagid[i] == 1){
    //     fprintf(screen,"\nx_COM[i] = %16.16f \nx_COM_use = %16.16f \nx[i][0] = %16.16f \ndelx = %16.16f",x_COM[i],x_COM_use,x[i][0]);
    //     fprintf(screen,"\ny_COM[i] = %16.16f \ny_COM_use = %16.16f \nx[i][1] = %16.16f \ndely = %16.16f",y_COM[i],y_COM_use,x[i][1]);
    //     fprintf(screen,"\nz_COM[i] = %16.16f \nz_COM_use = %16.16f \nx[i][2] = %16.16f \ndelz = %16.16f",z_COM[i],z_COM_use,x[i][2]);
    //   }
    // }

    delx = x[i][0]-x_COM_use;
    dely = x[i][1]-y_COM_use;
    delz = x[i][2]-z_COM_use;

    // additional terms from other confining potentials

    // add_factor_x = 2*strength*(force_add_factor*x_COM_use - x[i][0]);
    // add_factor_y = 2*strength*(force_add_factor*y_COM_use - x[i][1]);
    // add_factor_z = 2*strength*(force_add_factor*z_COM_use - x[i][2]);

    // f[i][0] += -2*force_factor*strength*(delx);
    // f[i][1] += -2*force_factor*strength*(dely);
    // f[i][2] += -2*force_factor*strength*(delz);

    // fprintf(screen,"\nx[i][0] = %16.16f",x[i][0]);
    // fprintf(screen,"\nx_COM[i] = %16.16f",x_COM[i]);
    // fprintf(screen,"\nx[i][1] = %16.16f",x[i][1]);
    // fprintf(screen,"\ny_COM[i] = %16.16f",y_COM[i]);
    // fprintf(screen,"\nx[i][2] = %16.16f",x[i][2]);
    // fprintf(screen,"\nz_COM[i] = %16.16f",z_COM[i]);

    f[i][0] += -2*strength*(delx);
    f[i][1] += -2*strength*(dely);
    f[i][2] += -2*strength*(delz);

    if (eflag) eng_vdwl += strength*(delx*delx + dely*dely + delz*delz);
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairConfine::allocate()
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

void PairConfine::settings(int narg, char **arg)
{
  if (narg != 4) error->all(FLERR,"Illegal pair_style command");

  strength = force->numeric(FLERR,arg[0]);
  cut_global = force->numeric(FLERR,arg[1]);
  box_len = force->numeric(FLERR,arg[2]);
  N_epe = force->numeric(FLERR,arg[3]);
  // fprintf(screen,"\nN_epe = %f",N_epe);

  force_factor = (N_epe-1.0)/N_epe;
  force_add_factor = (N_epe*N_epe - N_epe +1)/(N_epe*N_epe);
  // fprintf(screen,"\nforce_factor = %8.8f",force_factor);

  half_box_len = box_len/2;

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

void PairConfine::coeff(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR,"Incorrect args for pair coefficients");
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
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairConfine::init_style()
{
  neighbor->request(this,instance_me);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairConfine::init_one(int i, int j)
{
  if (setflag[i][j] == 0)
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);

  return cut[i][j];
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairConfine::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) fwrite(&cut[i][j],sizeof(double),1,fp);
    }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairConfine::read_restart(FILE *fp)
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

void PairConfine::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairConfine::read_restart_settings(FILE *fp)
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

void *PairConfine::extract(const char *str, int &dim)
{
  if (strcmp(str,"cut_global") == 0) {
    dim = 2;
    return (void *) &cut;
  }
  return NULL;
}