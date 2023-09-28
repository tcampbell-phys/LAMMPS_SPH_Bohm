#include "fix_fermi_no_cons.h"
#include <cstring>
#include "group.h"
#include "modify.h"
#include "error.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include <cmath>
#include <math.h>
#include "comm.h"
#include "neighbor.h"
#include "irregular.h"
#include "fix_deform.h"
#include "compute.h"
#include "kspace.h"
#include "respa.h"
#include "domain.h"
#include "memory.h"

using namespace LAMMPS_NS;
using namespace FixConst;

FixFermiNoCons::FixFermiNoCons(LAMMPS *lmp, int narg, char **arg):
	Fix(lmp, narg, arg), id_temp(NULL)
{
  if (narg < 6) error->all(FLERR,"Illegal fix Fermi command");

  int n = strlen(id) + 6;
  id_temp = new char[n];
  strcpy(id_temp,id);
  strcat(id_temp,"_temp");

  char **newarg = new char*[3];
  newarg[0] = id_temp;
  newarg[1] = group->names[igroup];
  newarg[2] = (char *) "temp";

  modify->add_compute(3,newarg);
  delete [] newarg;
  tcomputeflag = 1;

  t_target	= force->numeric(FLERR,arg[3]);
  full_box_len = force->numeric(FLERR,arg[4]);
  alpha_mass = force->numeric(FLERR,arg[5]);
  species_mass = force->numeric(FLERR,arg[6]);
  exp_cutoff = force->numeric(FLERR,arg[7]);
  mu = force->numeric(FLERR,arg[8]);

  alpha_half = 0.;
  alpha_minus_quarter = 0.;
  alpha_plus_quarter = 0.;
  alpha_dot = 0.;

}	

int FixFermiNoCons::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  return mask;
}

void FixFermiNoCons::init()
{
	dt = update->dt;
	// // fprintf(screen,"init() dt = %16.16f \n",dt);
	dthalf = 0.5 * update->dt;
	dt4 = 0.25 * update->dt;

	// to do: call boltz correct to units not just electron value
	boltz = 3.16681e-6;

	int icompute = modify->find_compute(id_temp);
	if (icompute < 0){
    	error->all(FLERR,"Temperature ID for fix fermi_nvt does not exist");
	}

	temperature = modify->compute[icompute];
}

void FixFermiNoCons::setup(int /*vflag*/)
{
	double** v = atom->v;
	// error->message(FLERR,"Checking FixFermiNoCons::setup() is called.");

	v_minus_quarter = v;
	v_plus_quarter = v;
}

void FixFermiNoCons::initial_integrate(int /*vflag*/)
{
	//update velocities by half a timestep, and positions by a full timestep:
	double** v = atom->v;
	double** x = atom->x;
	double** f = atom->f;

	int *mask = atom->mask;
	double potential_force;
	double friction_force;

	int nlocal = atom->nlocal;
	if (igroup == atom->firstgroup) nlocal = atom->nfirst;

	double beta;

	double force_multiplier = 0.;
	double friction_multiplier = 0.;
	double particle_energy = 0.;
	double exponent = 0.;
	double one_plus_exp = 0.;

	alpha_dot = 0.;

	double t_current;

	double ndim = 3;

	t_current = temperature->compute_scalar();
	//// fprintf(screen,"alpha_plus_quarter = %16.16f \n",alpha_plus_quarter);

	beta = 1 / (boltz * t_target);

	for(int i = 0; i < nlocal; ++i){
		if (mask[i] & groupbit){
			// // fprintf(screen,"i = %d \n",i);
			particle_energy = 0.;
			for(int j = 0; j < ndim; ++j){
				particle_energy += (v[i][j] * v[i][j]);
			// // fprintf(screen,"v[i][j] = %16.16f \n",v[i][j]);
			// // fprintf(screen,"in loop particle_energy = %16.16f \n",particle_energy);
			}
			particle_energy *= (0.5 * species_mass);
			//// fprintf(screen,"species_mass = %16.16f \n",species_mass);
			// // fprintf(screen,"post loop particle_energy = %16.16f \n",particle_energy);
			exponent = beta * (particle_energy - mu);
			// // fprintf(screen,"initial_integrate() exponent = %16.16f \n",exponent);
			one_plus_exp = 1. + exp(exponent);
			// // fprintf(screen,"initial_integrate() one_plus_exp = %16.16f \n",one_plus_exp);

			if(exponent > exp_cutoff){
				alpha_dot += 0.5 * alpha_mass * species_mass * (2. * particle_energy * beta	-3.);

			} 
			else {
				alpha_dot -= 0.5 * alpha_mass * species_mass * (2. * particle_energy * beta * 
					(2. - one_plus_exp) / one_plus_exp + 3.) * (one_plus_exp - 1.) / one_plus_exp;
			}
		}
		// fprintf(screen,"init alpha_dot = %16.16f \n",alpha_dot);
  }
	alpha_minus_quarter = alpha_plus_quarter;
	alpha_plus_quarter = alpha_minus_quarter + alpha_dot * dthalf;
	alpha_half = 0.5 * (alpha_minus_quarter + alpha_plus_quarter);

	//// fprintf(screen,"initial_integrate() alpha_half = %16.16f \n",alpha_half);

  	for(int i = 0; i < nlocal; ++i){
  		if (mask[i]&groupbit){
				particle_energy = 0.;
				for(int j = 0; j < ndim; ++j){
					particle_energy += (v[i][j] * v[i][j]);
				}
				particle_energy *= 0.5 * species_mass;
				exponent = beta * (particle_energy - mu);
				one_plus_exp = 1. + exp(exponent);

				if(exponent > exp_cutoff){
					force_multiplier = 1.;
					friction_multiplier = 0.5 * species_mass * beta;
				
				}
				else{
					force_multiplier = one_plus_exp * (log(one_plus_exp) - exponent);
					friction_multiplier = 0.5 * species_mass * beta * (one_plus_exp - 1.) / one_plus_exp;
				}
				for(int j = 0; j < ndim; ++j){
					v_minus_quarter[i][j] = v_plus_quarter[i][j];
					potential_force = f[i][j] * force_multiplier;
					// fprintf(screen,"\ninit f[i][j] = %16.16f",f[i][j]);
					friction_force = - friction_multiplier * alpha_half * v[i][j] * species_mass;
					v_plus_quarter[i][j] = v_minus_quarter[i][j] + dthalf * (potential_force + friction_force) / species_mass;
					v[i][j] = 0.5 * (v_minus_quarter[i][j] + v_plus_quarter[i][j]);
					x[i][j] += dt * v_plus_quarter[i][j];
			    // // fprintf(screen,"x[i][j] pre mod = %16.16f \n",x[i][j]);
		    	// x[i][j] = fmod(x[i][j],full_box_len);

			    // // fprintf(screen,"x[i][j] post mod = %16.16f \n",x[i][j]);
			}
    	}
    }

}

void FixFermiNoCons::final_integrate()
{
	// update velocities by another half timestep:
	double** v = atom->v;
	double** x = atom->x;
	double** f = atom->f;

	int *mask = atom->mask;
	double potential_force;
    double friction_force;

	int nlocal = atom->nlocal;
	if (igroup == atom->firstgroup) nlocal = atom->nfirst;

	double beta;

	double force_multiplier = 0.;
	double friction_multiplier = 0.;
	double particle_energy = 0.;
	double exponent = 0.;
	double exp_cutoff = 20.;
	double one_plus_exp = 0.;
	double t_current;

  alpha_dot = 0.;

	double ndim = 3;

	t_current = temperature->compute_scalar();
	//// fprintf(screen,"t_current = %16.16f \n",t_current);
	beta = 1 / (boltz * t_target);

	for(int i = 0; i < nlocal; ++i){
		if (mask[i] & groupbit){
			particle_energy = 0.;

			for(int j = 0; j < ndim; ++j){
				particle_energy += (v[i][j] * v[i][j]);
				//// fprintf(screen,"final_integrate() v[i][j] = %d \n",v[i][j]);
			}
			particle_energy *= 0.5 * species_mass;
			//// fprintf(screen,"initial_integrate() particle_energy = %16.16f \n",particle_energy);
			exponent = beta * (particle_energy - mu);
			//// fprintf(screen,"initial_integrate() exponent = %16.16f \n",exponent);
			one_plus_exp = 1. + exp(exponent);
			//// fprintf(screen,"initial_integrate() one_plus_exp = %16.16f \n",one_plus_exp);

			if(exponent > exp_cutoff){
				alpha_dot += 0.5 * alpha_mass * species_mass * (2. * particle_energy * beta	-3.);

			} 
			else {
				alpha_dot -= 0.5 * alpha_mass * species_mass * (2. * particle_energy * beta * 
					(2. - one_plus_exp) / one_plus_exp + 3.) * (one_plus_exp - 1.) / one_plus_exp;
			}
		}
		// fprintf(screen,"\nfinal alpha_dot = %16.16f \n",alpha_dot);
  }

	alpha_minus_quarter = alpha_plus_quarter;
	alpha_plus_quarter = alpha_minus_quarter + alpha_dot * dthalf;
	alpha_half = 0.5 * (alpha_minus_quarter + alpha_plus_quarter);

	//// fprintf(screen,"final_integrate() alpha_half = %16.16f \n",alpha_half);

  	for(int i = 0; i < nlocal; ++i){
  		if (mask[i] & groupbit){
				particle_energy = 0.;
				for(int j = 0; j < ndim; ++j){
					particle_energy += (v[i][j] * v[i][j]);
				}
				particle_energy *= 0.5 * species_mass;
				exponent = beta * (particle_energy - mu);
				// // fprintf(screen,"\nfinal exponent = %16.16f",exponent);
				one_plus_exp = 1. + exp(exponent);
				// // fprintf(screen,"\nfinal one_plus_exp = %16.16f",one_plus_exp);
				

				if(exponent > exp_cutoff){
					force_multiplier = 1.;
					friction_multiplier = 0.5 * species_mass * beta;
				}
				else{
					force_multiplier = one_plus_exp * (log(one_plus_exp) - exponent);
					friction_multiplier = 0.5 * species_mass * beta * (one_plus_exp - 1.) / one_plus_exp;
				}
				// fprintf(screen,"\nfinal force_multiplier = %16.16f",force_multiplier);
				
				for(int j = 0; j < ndim; ++j){

					v_minus_quarter[i][j] = v_plus_quarter[i][j];
		      potential_force = f[i][j] * force_multiplier;
					// fprintf(screen,"\nfinal f[i][j] = %16.16f",f[i][j]);
					// // fprintf(screen,"\nfinal potential_force = %16.16f",potential_force);
		      friction_force = - friction_multiplier * alpha_half * v[i][j] * species_mass;
		      // // fprintf(screen,"\nfinal friction_force = %16.16f",friction_force);
		      v_plus_quarter[i][j] = v_minus_quarter[i][j] + dthalf * (potential_force + friction_force) / species_mass;
					v[i][j] = 0.5 * (v_minus_quarter[i][j] + v_plus_quarter[i][j]);
					// // fprintf(screen,"\nfinal v[i][j] = %16.16f",v[i][j]);
	    	}
			}
  	}	
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write - for fix_fermi the only variable
   that needs to be transmitted forward is alpha_plus_quarter.
------------------------------------------------------------------------- */

void FixFermiNoCons::write_restart(FILE *fp)
{
  int nsize = 1;

  double *list;
  memory->create(list,nsize,"fermi:list");

  pack_restart_data(list);

  if (comm->me == 0) {
    int size = nsize * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(double),nsize,fp);
  }

  memory->destroy(list);
}

/* ----------------------------------------------------------------------
   pack restart data
------------------------------------------------------------------------- */

int FixFermiNoCons::pack_restart_data(double *list)
{
  int n = 0;
  list[n++] = alpha_plus_quarter;

  return n;
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixFermiNoCons::restart(char *buf)
{
  int n = 0;
  double *list = (double *) buf;
  alpha_plus_quarter = list[n++];
}
