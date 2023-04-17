#include "fermi_fix.h"
#include <cstring>
#include "group.h"
#include "modify.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

FixFermi::FixFermi(LAMMPS *lmp, int narg, char **arg):
	Fix(lmp, narg, arg)
{
  if (narg < 5)
  	  error->all(FLERR,"Illegal fix Fermi command");

  int n = strlen(id) + 6;
  id_temp = new char[n];
  strcpy(id_temp,id);
  strcat(id_temp,"_temp");

  char **newarg = new char*[3];
  newarg[0] = id_temp;
  newarg[1] = group->names[igroup];
  newarg[2] = (char *) "temp";
}	

int FixFermi::setmask()
{
  int mask = 0;
  mask |= FixConst::INITIAL_INTEGRATE;
  mask |= FixConst::FINAL_INTEGRATE;
  return mask;
}

void FixFermi::setup()
{
	double** v = atom->v;
	dt = update->dt;
	dthalf = 0.5 * update->dt;
	dt4 = 0.25 * update->dt;

	alpha_mass = 1.;
	alpha_half = 0.;
	alpha_minus_quarter = 0.;
	alpha_plus_quarter = 0.;
	alpha_dot = 0.;

	mu = 0.6622;

	v_minus_quarter = v;
	v_plus_quarter = v;

}

void FixFermi::initial_integrate()
{
	//update velocities by half a timestep, and positions by a full timestep (same procedure as in FixNH):
	double** v = atom->v;
	double** x = atom->x;
	double** f = atom->f;
	double* mass = atom->mass;
	double full_box_len = 26.4665;

	int nlocal = atom->nlocal;

	double boltz = force->boltz;
	double t_target = 40614.;
	double beta = 1 / (boltz * t_target);

	double force_multiplier = 0.;
	double friction_multiplier = 0.;
	double particle_energy = 0.;
	double exponent;
	double exp_cutoff = 20.;
	double one_plus_exp;

	double ndim = 3;



	for(int i = 0; i < nlocal; ++i){
		// need to check how to call species index of each particle:
		if(igroup = groupbit){
			continue;
		}
		for(int j = 0; j < ndim; ++j){
			particle_energy *= v[i][j] * v[i][j];
		}
		particle_energy *= 0.5 * mass[i];
		exponent = beta * (particle_energy - mu);
		one_plus_exp = 1. + exp(exponent);

		if(exponent > exp_cutoff){
			alpha_dot += 0.5 * alpha_mass * mass[i] * (2. * particle_energy * beta	-3.);

		} else {
			alpha_dot -= 0.5 * alpha_mass * mass[i] * (2. * particle_energy * beta * 
				(2. - one_plus_exp) / one_plus_exp + 3.) * (one_plus_exp - 1.) / one_plus_exp;
		}

		alpha_minus_quarter = alpha_plus_quarter;
		alpha_plus_quarter = alpha_minus_quarter + alpha_dot * dthalf;
		alpha_half = 0.5 * (alpha_minus_quarter + alpha_plus_quarter);
  }


  for(int i = 0; i < nlocal; ++i){
  	// need to check how to call species index of each particle:
  	if(igroup = groupbit){
			continue;
		}
		for(int j = 0; j < ndim; ++j){
			particle_energy *= v[i][j] * v[i][j];
		}
		particle_energy *= 0.5 * mass[i];
		exponent = beta * (particle_energy - mu);
		one_plus_exp = 1. + exp(exponent);

		if(exponent > exp_cutoff){
			force_multiplier = 1.;
			friction_multiplier = 0.5 * mass[i] * beta;
		
		}else{
			force_multiplier = one_plus_exp * (log(one_plus_exp) - exponent);
			friction_multiplier = 0.5 * mass[i] * beta * (one_plus_exp - 1.) / one_plus_exp;
		}
		for(int j = 0; j < ndim; ++j){
			v_minus_quarter[i][j] = v_plus_quarter[i][j];
      potential_force = forces[i][j] * force_multiplier;
      friction_force = - friction_multiplier * alpha_half * v[i][j] * mass[i];
      v_plus_quarter[i][j] = v_minus_quarter[i][j] + dthalf * (potential_force + 
      	friction_force) / mass[i];

      v[i][j] = 0.5 * (v_minus_quarter[i][j] + v_plus_quarter[i][j]);
      x[i][j] += dt * v_plus_quarter[i][j];
      x[i][j] = x[i][j] % full_box_len;
    }
  }
}

void FixFermi::final_integrate()
{
	// update velocities by another half timestep:
	double** v = atom->v;
	double** x = atom->x;
	double** f = atom->f;
	double* mass = atom->mass;
	double full_box_len = 33.2187;

	int nlocal = atom->nlocal;

	double boltz = force->boltz;
	double t_target = 40614.;
	double beta = 1 / (boltz * t_target);

	double force_multiplier = 0.;
	double friction_multiplier = 0.;
	double particle_energy = 0;
	double exponent;
	double exp_cutoff = 20.;
	double one_plus_exp;

	double ndim = 3;



	for(int i = 0; i < nlocal; ++i){
		// need to check how to call species index of each particle:
		if(igroup = groupbit){
			continue;
		}
		for(int j = 0; j < ndim; ++j){
			particle_energy *= v[i][j] * v[i][j];
		}
		particle_energy *= 0.5 * mass[i];
		exponent = beta * (particle_energy - mu);
		one_plus_exp = 1. + exp(exponent);

		if(exponent > exp_cutoff){
			alpha_dot += 0.5 * alpha_mass * mass[i] * (2. * particle_energy * beta	-3.);

		} else {
			alpha_dot -= 0.5 * alpha_mass * mass[i] * (2. * particle_energy * beta * 
				(2. - one_plus_exp) / one_plus_exp + 3.) * (one_plus_exp - 1.) / one_plus_exp;
		}

		alpha_minus_quarter = alpha_plus_quarter;
		alpha_plus_quarter = alpha_minus_quarter + alpha_dot * dthalf;
		alpha_half = 0.5 * (alpha_minus_quarter + alpha_plus_quarter);
  }


  for(int i = 0; i < nlocal; ++i){
  	// need to check how to call species index of each particle:
  	if(igroup = groupbit){
			continue;
		}
		for(int j = 0; j < ndim; ++j){
			particle_energy *= v[i][j] * v[i][j];
		}
		particle_energy *= 0.5 * mass[i];
		exponent = beta * (particle_energy - mu);
		one_plus_exp = 1. + exp(exponent);

		if(exponent > exp_cutoff){
			force_multiplier = 1.;
			friction_multiplier = 0.5 * mass[i] * beta;
		
		}else{
			force_multiplier = one_plus_exp * (log(one_plus_exp) - exponent);
			friction_multiplier = 0.5 * mass[i] * beta * (one_plus_exp - 1.) / one_plus_exp;
		}
		for(int j = 0; j < ndim; ++j){
			v_minus_quarter[i][j] = v_plus_quarter[i][j];
      potential_force = forces[i][j] * force_multiplier;
      friction_force = - friction_multiplier * alpha_half * v[i][j] * mass[i];
      v_plus_quarter[i][j] = v_minus_quarter[i][j] + dthalf * (potential_force + 
      	friction_force) / mass[i];

      v[i][j] = 0.5 * (v_minus_quarter[i][j] + v_plus_quarter[i][j]);
    }
  }
}