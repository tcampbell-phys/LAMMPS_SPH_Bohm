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

ComputeTempCentreMass::ComputeTempCentreMass(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg), CoM_vel(NULL), CoM_vel_all(NULL)
{
  if (narg != 3) error->all(FLERR,"Illegal compute temp_CoM command");

  // // fprintf(screen,"\nComputeTempCentreMass constructor...");

  CoM_vel = CoM_vel_all = NULL;

  scalar_flag = vector_flag = 1;
  size_vector = 6;
  extscalar = 0;
  extvector = 1;
  tempflag = 1;

  unallocated = 1;
  allocated = 0;
  vector = new double[size_vector];
 
}

/* ---------------------------------------------------------------------- */

ComputeTempCentreMass::~ComputeTempCentreMass()
{
  // fprintf(screen,"\nComputeTempCentreMass destructor...\n");
  if (!copymode)
  delete [] vector;
  if (allocated){
    deallocate();
  }
}

/* ---------------------------------------------------------------------- */

void ComputeTempCentreMass::setup()
{
  dynamic = 0;
  if (dynamic_user || group->dynamic[igroup]) dynamic = 1;
  dof_compute();
}

/* ---------------------------------------------------------------------- */

void ComputeTempCentreMass::dof_compute()
{
  adjust_dof_fix();
  natoms_temp = group->count(igroup);
  dof = domain->dimension * natoms_temp;
  dof -= extra_dof + fix_dof;
  if (dof > 0.0) tfactor = force->mvv2e / (dof * force->boltz);
  else tfactor = 0.0;
}

/* ---------------------------------------------------------------------- */

double ComputeTempCentreMass::compute_scalar()
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
  // // fprintf(screen,"\nComputeTempCentreMass::compute_scalar() tfactor = %f",tfactor);
  return scalar;
}

/* ---------------------------------------------------------------------- */

double ComputeTempCentreMass::compute_scalar_CoM(double N_epe, double particle_mass, int N_ele, int tag_ele_start)
{
  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int *tagid = atom->tag;

  double CoM_KE = 0.0;

  N_ele_compute_temp = N_ele;
  if (unallocated){
    deallocate();
    // fprintf(screen,"\nComputeTempCentreMass allocating CoM_vel memory...\n");
    // fprintf(screen,"\nComputeTempCentreMass particle_mass = %f\n",particle_mass);
    allocate();
    allocated = 1;
    unallocated = 0;
  }
  double tfactor_CoM;

  int ele_ind;

  int newton_pair = force->newton_pair;
  int size_CoM_vel = 3*N_ele_compute_temp;

  for (int i = 0; i < size_CoM_vel; i++) CoM_vel[i] = 0.0;

  //TODO: Calculate CoM Velocities and Communicate
  //careful not to double count contributions from electron CoMs

  for (int i = 0; i < nlocal; i++){
    if (mask[i] & groupbit){
      ele_ind = floor((tagid[i]-tag_ele_start)/N_epe);

      // if(ele_ind == 0){
      //   // fprintf(screen,"\ntagid[i] = %d",tagid[i]);
      //   // fprintf(screen,"\nv[i][0] = %f",v[i][0]);
      // }

      //vx 
      CoM_vel[ele_ind] += v[i][0]/N_epe;
      //vy
      CoM_vel[ele_ind+N_ele] += v[i][1]/N_epe;
      //vz
      CoM_vel[ele_ind+2*N_ele] += v[i][2]/N_epe;
    }
  }
  // // fprintf(screen,"\nA CoM_vel[0] = %f",CoM_vel[0]);
  MPI_Allreduce(CoM_vel,CoM_vel_all,size_CoM_vel,MPI_DOUBLE,MPI_SUM,world);
  // // fprintf(screen,"\nB CoM_vel_all[0] = %f",CoM_vel_all[0]);

  for (int i = 0; i < N_ele_compute_temp; i++){
    CoM_KE += (CoM_vel_all[i]*CoM_vel_all[i] + CoM_vel_all[i+N_ele]*CoM_vel_all[i+N_ele] + CoM_vel_all[i+2*N_ele]*CoM_vel_all[i+2*N_ele]) * particle_mass;
  }

  if (dynamic) dof_compute();
  if (dof < 0.0 && natoms_temp > 0.0)
    error->all(FLERR,"Temperature compute degrees of freedom < 0");
  tfactor_CoM = force->mvv2e / (((dof+extra_dof+fix_dof)/N_epe - (extra_dof+fix_dof)) * force->boltz);
  // fprintf(screen,"\nComputeTempCentreMass::compute_scalar_CoM() dof = %f",dof);
  // // fprintf(screen,"\nComputeTempCentreMass::compute_scalar_CoM() tfactor_CoM = %f",tfactor_CoM);
  CoM_KE *= tfactor_CoM;
  // fprintf(screen,"\nComputeTempCentreMass::compute_scalar_CoM() t_current = %f",CoM_KE);
  
  return CoM_KE;
}

/* ---------------------------------------------------------------------- */

void ComputeTempCentreMass::compute_vector()
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

/* ----------------------------------------------------------------------
   allocate memory that depends on # of electrons
------------------------------------------------------------------------- */

void ComputeTempCentreMass::allocate()
{
  // fprintf(screen,"\nComputeTempCentreMass::allocate()\n");
  // Centre of Mass array vx 0:N_ele, vy N_ele:2*N_ele, vz 2*N_ele:3_Nele
  CoM_vel = new double[3*N_ele_compute_temp]; 
  CoM_vel_all = new double[3*N_ele_compute_temp];
}

/* ----------------------------------------------------------------------
   deallocate memory that depends on # of electrons
------------------------------------------------------------------------- */

void ComputeTempCentreMass::deallocate()
{
  // fprintf(screen,"\nComputeTempCentreMass::deallocate()\n");
  delete [] CoM_vel;
  delete [] CoM_vel_all;
}