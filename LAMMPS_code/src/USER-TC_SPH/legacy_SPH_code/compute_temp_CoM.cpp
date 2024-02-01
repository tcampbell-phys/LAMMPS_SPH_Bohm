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

#include "compute_temp_CoM.h"
#include <mpi.h>
#include "atom.h"
#include "update.h"
#include "force.h"
#include "domain.h"
#include "group.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeTempCoM::ComputeTempCoM(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 3) error->all(FLERR,"Illegal compute temp command");

  scalar_flag = vector_flag = 1;
  // size_vector = 6;
  // size_vector_CoM_vel = 1;
  extscalar = 0;
  extvector = 1;
  tempflag = 1;
  nmax = 0;

  if (narg == 4){
    fprintf(screen,"\nin CoM...");
    N_ele_compute_temp = force->numeric(FLERR,arg[3]);
    size_vector = 3*N_ele_compute_temp;
    fprintf(screen,"\nComputeTempCoM size_vector = %d",size_vector);
    CoM_vel_flag = 1;
    // vector_CoM_vel = new double[size_vector_CoM_vel];
  } else if (narg ==3){
    fprintf(screen,"\nin standard...");
    size_vector = 6;
    fprintf(screen,"\nComputeTempCoM size_vector = %d",size_vector);
  }

  fprintf(screen,"\nComputeTempCoM constructor...");
  // printf("\n%s",arg[2]);

  // if (narg == 4){
  //   fprintf(screen,"\nin CoM...");
  //   N_ele_compute_temp = force->numeric(FLERR,arg[3]);
  //   size_vector_CoM_vel = 3*N_ele_compute_temp;
  //   fprintf(screen,"\nComputeTempCoM size_vector_CoM_vel = %d",size_vector_CoM_vel);
  //   CoM_vel_flag = 1;
  //   vector_CoM_vel = new double[size_vector_CoM_vel];
  // } else if (narg ==3){
  //   fprintf(screen,"\nin standard...");
  //   fprintf(screen,"\nComputeTempCoM size_vector_CoM_vel = %d",size_vector_CoM_vel);
  // }

  vector = new double[size_vector];
 
}

/* ---------------------------------------------------------------------- */

ComputeTempCoM::~ComputeTempCoM()
{
  fprintf(screen,"\nComputeTempCoM destructor...");
  if (!copymode)
    // fprintf(screen,"\nCoM_vel_flag = %d",CoM_vel_flag);
    delete [] vector;
    // if (CoM_vel_flag){
    //   fprintf(screen,"\nCoM_vel_flag = 1");
    //   delete [] vector_CoM_vel;
    // }
}

/* ---------------------------------------------------------------------- */

void ComputeTempCoM::setup()
{
  dynamic = 0;
  if (dynamic_user || group->dynamic[igroup]) dynamic = 1;
  dof_compute();
}

/* ---------------------------------------------------------------------- */

void ComputeTempCoM::dof_compute()
{
  adjust_dof_fix();
  natoms_temp = group->count(igroup);
  dof = domain->dimension * natoms_temp;
  dof -= extra_dof + fix_dof;
  if (dof > 0.0) tfactor = force->mvv2e / (dof * force->boltz);
  else tfactor = 0.0;
}

/* ---------------------------------------------------------------------- */

double ComputeTempCoM::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double t = 0.0;

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
        t += (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]) * rmass[i];
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
        t += (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]) *
          mass[type[i]];
  }

  MPI_Allreduce(&t,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  if (dynamic) dof_compute();
  if (dof < 0.0 && natoms_temp > 0.0)
    error->all(FLERR,"Temperature compute degrees of freedom < 0");
  scalar *= tfactor;
  // fprintf(screen,"\nComputeTempCoM::compute_scalar() tfactor = %f",tfactor);
  return scalar;
}

/* ---------------------------------------------------------------------- */

double ComputeTempCoM::compute_scalar_CoM(double N_epe, int N_ele, int tag_ele_start)
{
  invoked_scalar = update->ntimestep;
  invoked_vector = update->ntimestep;

  // int vector_size = 3*N_ele;

  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int *tagid = atom->tag;

  double t = 0.0;
  double v_CoM[size_vector]; // Centre of Mass array vx 0:N_ele, vy N_ele:2*N_ele, vz 2*N_ele:3_Nele
  double tfactor_CoM;

  int ele_ind;

  int newton_pair = force->newton_pair;

  for (int i = 0; i < size_vector; i++) v_CoM[i] = 0.0;

  //TODO: Calculate CoM Velocities and Communicate
  //careful not to double count contributions from electron CoMs

  for (int i = 0; i < nlocal; i++){
    if (mask[i] & groupbit){
      // fprintf(screen,"\ntagid[i] = %d",tagid[i]);
      ele_ind = floor((tagid[i]-tag_ele_start)/N_epe);

      // fprintf(screen,"\nele_ind = %d",ele_ind);
      //vx 
      v_CoM[ele_ind] += v[i][0]/N_epe;
      //vy
      v_CoM[ele_ind+N_ele] += v[i][1]/N_epe;
      //vz
      v_CoM[ele_ind+2*N_ele] += v[i][2]/N_epe;
    }
  }
  fprintf(screen,"\nA v_CoM[0] = %f",v_CoM[0]);
  MPI_Allreduce(v_CoM,vector,size_vector,MPI_DOUBLE,MPI_SUM,world);
  fprintf(screen,"\nB vector[0] = %f",vector[0]);

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
        t += (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]) * rmass[i];
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
        t += (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]) *
          mass[type[i]];
  }

  MPI_Allreduce(&t,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  if (dynamic) dof_compute();
  if (dof < 0.0 && natoms_temp > 0.0)
    error->all(FLERR,"Temperature compute degrees of freedom < 0");
  tfactor_CoM = force->mvv2e / (((dof+extra_dof+fix_dof)/N_epe - (extra_dof+fix_dof)) * force->boltz);
  // fprintf(screen,"\nComputeTempCoM::compute_scalar_CoM() tfactor_CoM = %f",tfactor_CoM);
  scalar *= tfactor_CoM;
  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeTempCoM::compute_vector()
{
  int i;

  invoked_vector = update->ntimestep;

  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double massone,t[6];
  for (i = 0; i < 6; i++) t[i] = 0.0;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      t[0] += massone * v[i][0]*v[i][0];
      t[1] += massone * v[i][1]*v[i][1];
      t[2] += massone * v[i][2]*v[i][2];
      t[3] += massone * v[i][0]*v[i][1];
      t[4] += massone * v[i][0]*v[i][2];
      t[5] += massone * v[i][1]*v[i][2];
    }

  MPI_Allreduce(t,vector,6,MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < 6; i++) vector[i] *= force->mvv2e;
}

// /* ---------------------------------------------------------------------- */

// int ComputeTempCoM::pack_reverse_comm(int n, int first, double *buf)
// {
//   int i,m,last;

//   m = 0;
//   last = first + n;

//   for (i = first; i < last; i++){
//     buf[m++] = vx_CoM[i];
//     buf[m++] = vy_CoM[i];
//     buf[m++] = vz_CoM[i];
//   }
//   return m;
// }

// /* ---------------------------------------------------------------------- */

// void ComputeTempCoM::unpack_reverse_comm(int n, int *list, double *buf)
// {
//   int i,j,m;

//   m = 0;

//   for (i = 0; i < n; i++) {
//     j = list[i];
//     vx_CoM[j] += buf[m++];
//     vy_CoM[j] += buf[m++];
//     vz_CoM[j] += buf[m++];
//   }
// }
