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
  fprintf(screen,"\nconfine constructor...");
}

/* ---------------------------------------------------------------------- */

PairConfine::~PairConfine()
{
  fprintf(screen,"\nconfine destructor...");
}

/* ---------------------------------------------------------------------- */

void PairConfine::compute(int eflag, int vflag)
{
  fprintf(screen,"\nconfine compute...");
  int i,ii,inum;
  double e_confine;
  double delx,dely,delz;
  int *ilist;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  fprintf(screen,"A");
  double *x_COM = atom->x_COM;
  double *y_COM = atom->y_COM;
  double *z_COM = atom->z_COM;
  fprintf(screen,"B");
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];

    delx = x[i][0]-x_COM[i];
    dely = x[i][1]-y_COM[i];
    delz = x[i][2]-z_COM[i];

    f[i][0] += -2*strength*(delx);
    f[i][1] += -2*strength*(dely);
    f[i][2] += -2*strength*(delz);

    if (eflag) e_confine = strength*(delx*delx + dely*dely + delz*delz);
    if (evflag) ev_tally(i,0,nlocal,1,
                          e_confine,0.0,0.0,0.0,0.0,0.0);

  }
}
/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairConfine::allocate()
{
  allocated = 1;
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairConfine::settings(int narg, char **arg)
{
  fprintf(screen,"\nconfine settings A...");
  if (narg != 1) error->all(FLERR,"Illegal pair_style command");
  fprintf(screen,"\nconfine settings B...");
  strength = force->numeric(FLERR,arg[0]);
  fprintf(screen,"\nconfine settings C...");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairConfine::coeff(int narg, char **arg)
{
  fprintf(screen,"\nconfine coeff...");
  if (narg != 2) error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairConfine::init_style()
{
  fprintf(screen,"\nconfine init_style...");
  neighbor->request(this,instance_me);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairConfine::init_one(int i, int j)
{
  fprintf(screen,"\nconfine init_one...");
  return strength;
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairConfine::write_restart(FILE *fp)
{
  fprintf(screen,"\nconfine write_restart...");
  write_restart_settings(fp);

  fwrite(&strength,sizeof(int),1,fp);
      
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairConfine::read_restart(FILE *fp)
{
  fprintf(screen,"\nconfine read_restart...");
  
  read_restart_settings(fp);

  int me = comm->me;
  if (me == 0) utils::sfread(FLERR,&strength,sizeof(int),1,fp,NULL,error);
  MPI_Bcast(&strength,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairConfine::write_restart_settings(FILE *fp)
{
  fprintf(screen,"\nconfine write_restart_settings...");
  
  fwrite(&strength,sizeof(double),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairConfine::read_restart_settings(FILE *fp)
{
  fprintf(screen,"\nconfine read_restart_settings...");
  
  if (comm->me == 0) {
    utils::sfread(FLERR,&strength,sizeof(double),1,fp,NULL,error);
  }
  MPI_Bcast(&strength,1,MPI_DOUBLE,0,world);
}

/* ---------------------------------------------------------------------- */

void *PairConfine::extract(const char *str, int &dim)
{
  fprintf(screen,"\nconfine extract...");
  
  return NULL;
}