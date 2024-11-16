/* -*- c++ -*- ----------------------------------------------------------------------
Thomas Campbell (Oxford)
------------------------------------------------------------------------- */

#include "pair_confine_full.h"
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

PairConfineFull::PairConfineFull(LAMMPS *lmp) : Pair(lmp)
{
  boltz_val = force->boltz;
  hbar_val = (force->hplanck)/(2*M_PI);

  all_delx = NULL;
  all_dely = NULL;
  all_delz = NULL;
  unallocated = 1;
  // // fprintf(screen,"\nIn pair_confine with correct forcing term...");
}

/* ---------------------------------------------------------------------- */

PairConfineFull::~PairConfineFull()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    deallocate();
  }
}

/* ---------------------------------------------------------------------- */

void PairConfineFull::compute(int eflag, int vflag)
{
  // fprintf(screen,"\nCompute 0");
  int i,ii,inum,itype;
  double delx,dely,delz;
  int *ilist;

  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;

  double e_confine;

  int *tagid = atom->tag;

  // // fprintf(screen,"\nCompute 1");
  // int ntimestep = update->ntimestep;

  double *x_COM = atom->x_COM;
  double *y_COM = atom->y_COM;
  double *z_COM = atom->z_COM;

  // // fprintf(screen,"\nCompute 2");

  double x_COM_use,y_COM_use,z_COM_use;

  // // fprintf(screen,"\nCompute 3");

  int ele_ind;
  int tag_val;
  int tag_start;
  // // fprintf(screen,"\nCompute 4");

  // // fprintf(screen,"\nCompute A");

  // if (unallocated){
  //   deallocate();
  //   // // // fprintf(screen,"\nFixDynamicWidthsCoM allocating mu tau memory...\n");
  //   allocate_array();
  //   unallocated = 0;
  // }

  for (int a = 0; a < N_SPH; a++){ 
    ind_delx[a] = 0.0;
    ind_dely[a] = 0.0;
    ind_delz[a] = 0.0;
  }

  // // fprintf(screen,"\nCompute B");

  // // fprintf(screen,"\nCompute B2");

  inum = list->inum;
  ilist = list->ilist;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];

    // Choose correct projection of Centre of Mass

    if (x[i][0]-x_COM[i] > half_box_len){
      // // fprintf(screen,"\nBasic Ax");
      x_COM_use = x_COM[i] + box_len;
    } else if (x[i][0]-x_COM[i] < -half_box_len){
      // // fprintf(screen,"\nBasic Bx");
      x_COM_use = x_COM[i] - box_len;
    } else{
      x_COM_use = x_COM[i];
    }
    if (x[i][1]-y_COM[i] > half_box_len){
      // // fprintf(screen,"\nBasic Ay");
      y_COM_use = y_COM[i] + box_len;
    } else if (x[i][1]-y_COM[i] < -half_box_len){
      // // fprintf(screen,"\nBasic By");
      y_COM_use = y_COM[i] - box_len;
    } else{
      y_COM_use = y_COM[i];
    }
    if (x[i][2]-z_COM[i] > half_box_len){
      // // fprintf(screen,"\nBasic Az");
      z_COM_use = z_COM[i] + box_len;
    } else if (x[i][2]-z_COM[i] < -half_box_len){
      // // fprintf(screen,"\nBasic Bz");
      z_COM_use = z_COM[i] - box_len;
    } else{
      z_COM_use = z_COM[i];
    }

    delx = x[i][0]-x_COM_use;
    dely = x[i][1]-y_COM_use;
    delz = x[i][2]-z_COM_use;

    ind_delx[tagid[i]-tag_ele_start] = delx;
    ind_dely[tagid[i]-tag_ele_start] = dely;
    ind_delz[tagid[i]-tag_ele_start] = delz;

    // // fprintf(screen,"\ntagid[i] = %d",tagid[i]);
    // // fprintf(screen,"\nind_delx[tagid[i]-tag_ele_start] = %16.16f",ind_delx[tagid[i]-tag_ele_start]);
    // // fprintf(screen,"\nind_dely[tagid[i]-tag_ele_start] = %16.16f",ind_dely[tagid[i]-tag_ele_start]);
    // // fprintf(screen,"\nind_delz[tagid[i]-tag_ele_start] = %16.16f",ind_delz[tagid[i]-tag_ele_start]);
  }

  // communicate delx/y/z terms globally

  MPI_Allreduce(ind_delx,all_delx,N_SPH,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(ind_dely,all_dely,N_SPH,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(ind_delz,all_delz,N_SPH,MPI_DOUBLE,MPI_SUM,world);

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];

    // fprintf(screen,"\ntagid[i] = %d",tagid[i]);
    // fprintf(screen,"\ntagid[i]- tag_ele_start= %d",tagid[i]-tag_ele_start);
    // fprintf(screen,"\narray_ind = %d",tagid[i]-tag_ele_start);
    // // fprintf(screen,"\nall_delx[tagid[i]-tag_ele_start] = %16.16f",ind_delx[tagid[i]-tag_ele_start]);
    // // fprintf(screen,"\nall_dely[tagid[i]-tag_ele_start] = %16.16f",ind_dely[tagid[i]-tag_ele_start]);
    // // fprintf(screen,"\nall_delz[tagid[i]-tag_ele_start] = %16.16f",ind_delz[tagid[i]-tag_ele_start]);
  

    f[i][0] += -2*strength*force_factor*(all_delx[tagid[i]-tag_ele_start]);
    f[i][1] += -2*strength*force_factor*(all_dely[tagid[i]-tag_ele_start]);
    f[i][2] += -2*strength*force_factor*(all_delz[tagid[i]-tag_ele_start]);


    ele_ind = floor((tagid[i]-tag_ele_start)/N_epe);
    // fprintf(screen,"\nele_ind = %d",ele_ind);
    tag_start = (ele_ind*N_epe);
    // fprintf(screen,"\ntag_start = %d",tag_start);

    // additional terms

    for (int j = 0; j < N_epe;j++){
      tag_val = tag_start+j;
      // fprintf(screen,"\ntag_val = %d",tag_val);
      if (tag_val == (tagid[i]-tag_ele_start)){
        // fprintf(screen,"\ncontinuing...");
        continue;
      }
      f[i][0] += 2*strength*all_delx[tag_val]/N_epe;
      f[i][1] += 2*strength*all_dely[tag_val]/N_epe;
      f[i][2] += 2*strength*all_delz[tag_val]/N_epe;
    }

    if (eflag) eng_vdwl += strength*(all_delx[tagid[i]-tag_ele_start]*all_delx[tagid[i]-tag_ele_start] + all_dely[tagid[i]-tag_ele_start]*all_dely[tagid[i]-tag_ele_start] + all_delz[tagid[i]-tag_ele_start]*all_delz[tagid[i]-tag_ele_start]);
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairConfineFull::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(cut,n+1,n+1,"pair:cut");

  all_delx = new double[N_SPH]; 
  all_dely = new double[N_SPH];
  all_delz = new double[N_SPH];
  ind_delx = new double[N_SPH]; 
  ind_dely = new double[N_SPH];
  ind_delz = new double[N_SPH];
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairConfineFull::settings(int narg, char **arg)
{
  if (narg != 6) error->all(FLERR,"Illegal pair_style command");
  // fprintf(screen,"\nSettings A");
  strength = force->numeric(FLERR,arg[0]);
  cut_global = force->numeric(FLERR,arg[1]);
  box_len = force->numeric(FLERR,arg[2]);
  N_epe = force->numeric(FLERR,arg[3]);
  tag_ele_start = force->numeric(FLERR,arg[4]);
  N_SPH = force->numeric(FLERR,arg[5]);
  // // fprintf(screen,"\nN_epe = %f",N_epe);

  // fprintf(screen,"\nSettings B");

  force_factor = (N_epe-1)/N_epe;
  force_add_factor = (N_epe*N_epe - N_epe +1)/(N_epe*N_epe);
  // // fprintf(screen,"\nforce_factor = %8.8f",force_factor);

  // fprintf(screen,"\nSettings C");

  half_box_len = box_len/2;

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }

  // fprintf(screen,"\nSettings D");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairConfineFull::coeff(int narg, char **arg)
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

void PairConfineFull::init_style()
{
  neighbor->request(this,instance_me);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairConfineFull::init_one(int i, int j)
{
  if (setflag[i][j] == 0)
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);

  return cut[i][j];
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairConfineFull::write_restart(FILE *fp)
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

void PairConfineFull::read_restart(FILE *fp)
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

void PairConfineFull::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairConfineFull::read_restart_settings(FILE *fp)
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

void *PairConfineFull::extract(const char *str, int &dim)
{
  if (strcmp(str,"cut_global") == 0) {
    dim = 2;
    return (void *) &cut;
  }
  return NULL;
}

/* ----------------------------------------------------------------------
   deallocate memory that depends on # of electrons
------------------------------------------------------------------------- */

void PairConfineFull::deallocate()
{
  // // fprintf(screen,"\nFixDynamicWidthsCoM::deallocate()\n");
  delete [] all_delx;
  delete [] all_dely;
  delete [] all_delz;
  delete [] ind_delx;
  delete [] ind_dely;
  delete [] ind_delz;

}