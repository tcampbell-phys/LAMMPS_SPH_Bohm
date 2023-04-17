/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Edited by Thomas Campbell 2021.
------------------------------------------------------------------------- */

#include "fix_nve_cons.h"
#include <cstring>
#include "atom.h"
#include "force.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "memory.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixNVEcons::FixNVEcons(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (strcmp(style,"nve/sphere") != 0 && narg < 3)
    error->all(FLERR,"Illegal fix nve command");

  dynamic_group_allow = 1;
  time_integrate = 1;

  int ndim = 3;
  old_nlocal = atom->nlocal;

  memory->create(f_corrected,atom->nlocal,ndim,"fix:f_corrected");
}

/* ---------------------------------------------------------------------- */

FixNVEcons::~FixNVEcons()
{
memory->destroy(f_corrected);
}

/* ---------------------------------------------------------------------- */

int FixNVEcons::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= INITIAL_INTEGRATE_RESPA;
  mask |= FINAL_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixNVEcons::init()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;

  if (strstr(update->integrate_style,"respa"))
    step_respa = ((Respa *) update->integrate)->step;
  fprintf(screen,"Calling fix_nve_cons.cpp...\n");

}

/* ----------------------------------------------------------------------
   allow for both per-type and per-atom mass
------------------------------------------------------------------------- */

void FixNVEcons::initial_integrate(int /*vflag*/)
{
  double dtfm;

  // update v and x of atoms in group

  //to ensure energy conservation for Bohm GWP with de Broglie wavepacket widths,
  //impose beta correction to the force calculations

  double *beta_correction = atom->beta_correction;

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int ndim = 3;

  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  if (nlocal != old_nlocal){
		memory->grow(f_corrected,nlocal,ndim,"fix:f_corrected");
		old_nlocal = nlocal;
	}

	for(int i = 0; i < nlocal; ++i){
		if (mask[i] & groupbit){
			//fprintf(screen,"pre assign initial f_corrected [i,x] = %.9f \n",f_corrected[i][0]);
			//fprintf(screen,"pre assign inital f [i,x] = %.9f \n",f[i][0]);
			//fprintf(screen,"pre assign pos[i][x] = %.9f \n",x[i][0]);
			//fprintf(screen,"pre assign pos[i][y] = %.9f \n",x[i][1]);
			//fprintf(screen,"pre assign pos[i][z] = %.9f \n",x[i][2]);
			//fprintf(screen,"pre initial calc f_corrected [i,x] = %.9f \n",f_corrected[i][0]);
			//fprintf(screen,"pre initial calc f [i,x] = %.9f \n",f[0][0]);
			//fprintf(screen,"pre inital beta_correction [i] = %.15f \n",beta_correction[i]);
			for(int j = 0; j < ndim; ++j){
				f_corrected[i][j] = f[i][j];
			}
		}
	}

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / rmass[i];

        //fprintf(screen,"pre initial calc f [i,x] = %.15f \n",f[i][0]);
			  //fprintf(screen,"pre inital f_corrected [i,X] = %.15f \n",f_corrected[i][0]);

        f_corrected[i][0] *= 1/(1-beta_correction[i]);
        f_corrected[i][1] *= 1/(1-beta_correction[i]);
        f_corrected[i][2] *= 1/(1-beta_correction[i]);

        //fprintf(screen,"post initial calc f [i,x] = %.15f \n",f[i][0]);
			  //fprintf(screen,"post inital f_corrected [i,X] = %.15f \n",f_corrected[i][0]);

        v[i][0] += dtfm * f_corrected[i][0];
        v[i][1] += dtfm * f_corrected[i][1];
        v[i][2] += dtfm * f_corrected[i][2];
        x[i][0] += dtv * v[i][0];
        x[i][1] += dtv * v[i][1];
        x[i][2] += dtv * v[i][2];
      }

  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / mass[type[i]];

        //fprintf(screen,"pre initial calc f [i,x] = %.15f \n",f[i][0]);
			  //fprintf(screen,"pre inital f_corrected [i,X] = %.15f \n",f_corrected[i][0]);

        f_corrected[i][0] *= 1/(1-beta_correction[i]);
        f_corrected[i][1] *= 1/(1-beta_correction[i]);
        f_corrected[i][2] *= 1/(1-beta_correction[i]);

        //fprintf(screen,"post initial calc f [i,x] = %.15f \n",f[i][0]);
			  //fprintf(screen,"post inital f_corrected [i,X] = %.15f \n",f_corrected[i][0]);

        v[i][0] += dtfm * f_corrected[i][0];
        v[i][1] += dtfm * f_corrected[i][1];
        v[i][2] += dtfm * f_corrected[i][2];
        x[i][0] += dtv * v[i][0];
        x[i][1] += dtv * v[i][1];
        x[i][2] += dtv * v[i][2];
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixNVEcons::final_integrate()
{
  double dtfm;

  // update v of atoms in group

  //to ensure energy conservation for Bohm GWP with de Broglie wavepacket widths,
  //impose beta correction to the force calculations
  
  double *beta_correction = atom->beta_correction;

  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int ndim = 3;

  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  if (nlocal != old_nlocal){
		memory->grow(f_corrected,nlocal,ndim,"fix:f_corrected");
		old_nlocal = nlocal;
	}

	for(int i = 0; i < nlocal; ++i){
		if (mask[i] & groupbit){
			//fprintf(screen,"pre assign initial f_corrected [i,x] = %.9f \n",f_corrected[i][0]);
			//fprintf(screen,"pre assign inital f [i,x] = %.9f \n",f[i][0]);
			//fprintf(screen,"pre assign pos[i][x] = %.9f \n",x[i][0]);
			//fprintf(screen,"pre assign pos[i][y] = %.9f \n",x[i][1]);
			//fprintf(screen,"pre assign pos[i][z] = %.9f \n",x[i][2]);
			//fprintf(screen,"pre initial calc f_corrected [i,x] = %.9f \n",f_corrected[i][0]);
			//fprintf(screen,"pre initial calc f [i,x] = %.9f \n",f[0][0]);
			//fprintf(screen,"pre inital beta_correction [i] = %.15f \n",beta_correction[i]);
			for(int j = 0; j < ndim; ++j){
				f_corrected[i][j] = f[i][j];
			}
		}
	}

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / rmass[i];

        f_corrected[i][0] *= 1/(1-beta_correction[i]);
        f_corrected[i][1] *= 1/(1-beta_correction[i]);
        f_corrected[i][2] *= 1/(1-beta_correction[i]);

        v[i][0] += dtfm * f_corrected[i][0];
        v[i][1] += dtfm * f_corrected[i][1];
        v[i][2] += dtfm * f_corrected[i][2];
      }

  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / mass[type[i]];
        f_corrected[i][0] *= 1/(1-beta_correction[i]);
        f_corrected[i][1] *= 1/(1-beta_correction[i]);
        f_corrected[i][2] *= 1/(1-beta_correction[i]);

        v[i][0] += dtfm * f_corrected[i][0];
        v[i][1] += dtfm * f_corrected[i][1];
        v[i][2] += dtfm * f_corrected[i][2];
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixNVEcons::initial_integrate_respa(int vflag, int ilevel, int /*iloop*/)
{
  dtv = step_respa[ilevel];
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;

  // innermost level - NVE update of v and x
  // all other levels - NVE update of v

  if (ilevel == 0) initial_integrate(vflag);
  else final_integrate();
}

/* ---------------------------------------------------------------------- */

void FixNVEcons::final_integrate_respa(int ilevel, int /*iloop*/)
{
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;
  final_integrate();
}

/* ---------------------------------------------------------------------- */

void FixNVEcons::reset_dt()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
}
