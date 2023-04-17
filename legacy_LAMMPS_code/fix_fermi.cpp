#include "fix_fermi.h"
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

FixFermi::FixFermi(LAMMPS *lmp, int narg, char **arg):
	Fix(lmp, narg, arg), id_temp(NULL)
{
  if (narg < 7) error->all(FLERR,"Illegal fix Fermi command");

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
  initial_alpha = force->numeric(FLERR,arg[9]);

  alpha_half = 0.;
  alpha_minus_quarter = 0.;
  alpha_plus_quarter = initial_alpha;
  alpha_dot = 0.;

  int ndim = 3;
  old_nlocal = atom->nlocal;

  memory->create(f_corrected,atom->nlocal,ndim,"fix:f_corrected");


}	

FixFermi::~FixFermi()
{
memory->destroy(f_corrected);
}

int FixFermi::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  return mask;
}

void FixFermi::init()
{
	dt = update->dt;
	//fprintf(screen,"init() dt = %f \n",dt);
	dthalf = 0.5 * update->dt;
	dt4 = 0.25 * update->dt;

	boltz = 3.16681e-6;

	int icompute = modify->find_compute(id_temp);
	if (icompute < 0){
    	error->all(FLERR,"Temperature ID for fix fermi_nvt does not exist");
	}

	fprintf(screen,"alpha_plus_quarter = %f\n",alpha_plus_quarter);

	temperature = modify->compute[icompute];
}

void FixFermi::setup(int /*vflag*/)
{
	double** v = atom->v;
	//error->message(FLERR,"Checking FixFermi::setup() is called.");

    v_minus_quarter = v;
    v_plus_quarter = v;
}

void FixFermi::initial_integrate(int /*vflag*/)
{	
	//fprintf(screen,"alpha_plus_quarter = %f",alpha_plus_quarter);
	//update velocities by half a timestep, and positions by a full timestep (same procedure as in FixFermi):
	double** v = atom->v;
	double** x = atom->x;
	double** f = atom->f;

	int *mask = atom->mask;
	double potential_force;
	double friction_force;

	//to ensure energy conservation for Bohm GWP with de Broglie wavepacket widths,
	//impose beta correction to the force calculations

	double *beta_correction = atom->beta_correction;
	

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



	t_current = temperature->compute_scalar();
	//fprintf(screen,"alpha_plus_quarter = %f \n",alpha_plus_quarter);

	beta = 1 / (boltz * t_target);

	for(int i = 0; i < nlocal; ++i){
		if (mask[i] & groupbit){
			//fprintf(screen,"i = %d \n",i);
			particle_energy = 0.;
			for(int j = 0; j < ndim; ++j){
				particle_energy += (v[i][j] * v[i][j]);
			//	fprintf(screen,"v[i][j] = %f \n",v[i][j]);
			//	fprintf(screen,"in loop particle_energy = %f \n",particle_energy);
			}
			particle_energy *= (0.5 * species_mass);
			//fprintf(screen,"species_mass = %f \n",species_mass);
			//fprintf(screen,"post loop particle_energy = %f \n",particle_energy);
			exponent = beta * (particle_energy - mu);
			//fprintf(screen,"initial_integrate() exponent = %f \n",exponent);
			one_plus_exp = 1. + exp(exponent);
			//fprintf(screen,"initial_integrate() one_plus_exp = %f \n",one_plus_exp);

			if(exponent > exp_cutoff){
				alpha_dot += 0.5 * alpha_mass * species_mass * (2. * particle_energy * beta	-3.);

			} 
			else {
				alpha_dot -= 0.5 * alpha_mass * species_mass * (2. * particle_energy * beta * 
					(2. - one_plus_exp) / one_plus_exp + 3.) * (one_plus_exp - 1.) / one_plus_exp;
			}
			//fprintf(screen,"alpha_dot = %f \n",alpha_dot);
		}
    }
	alpha_minus_quarter = alpha_plus_quarter;
	alpha_plus_quarter = alpha_minus_quarter + alpha_dot * dthalf;
	alpha_half = 0.5 * (alpha_minus_quarter + alpha_plus_quarter);

	//fprintf(screen,"initial_integrate() alpha_half = %f \n",alpha_half);

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
				//apply beta_correction
				//fprintf(screen,"f_corrected[i][j] pre mod = %f \n",f_corrected[i][j]);
				//fprintf(screen,"f [i,j] pre mod = %.9f \n",f[i][j]);
				f_corrected[i][j] *= 1/(1-beta_correction[i]);
				//fprintf(screen,"f_corrected [i,j] post mod = %.9f \n",f_corrected[i][j]);
				//fprintf(screen,"f [i,j] post mod = %.9f \n",f[i][j]);
				//fprintf(screen,"f_corrected[i][j] post mod = %f \n",f_corrected[i][j]);

				v_minus_quarter[i][j] = v_plus_quarter[i][j];
	      		potential_force = f_corrected[i][j] * force_multiplier;
	      		friction_force = - friction_multiplier * alpha_half * v[i][j] * species_mass;
	     	 	v_plus_quarter[i][j] = v_minus_quarter[i][j] + dthalf * (potential_force + friction_force) / species_mass;
				v[i][j] = 0.5 * (v_minus_quarter[i][j] + v_plus_quarter[i][j]);
			    x[i][j] += dt * v_plus_quarter[i][j];
			    //fprintf(screen,"x[i][j] pre mod = %f \n",x[i][j]);
			    x[i][j] = fmod(x[i][j],full_box_len);

			    //fprintf(screen,"x[i][j] post mod = %f \n",x[i][j]);
			}
    	}
    }
	//fprintf(screen,"post initial calc f_corrected [0,x] = %.9f \n",f_corrected[0][0]);
	//fprintf(screen,"post initial calc f [0,x] = %.9f \n",f[0][0]);
}

void FixFermi::final_integrate()
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

	//to ensure energy conservation for Bohm GWP with de Broglie wavepacket widths,
	//impose beta correction to the force calculations

	double *beta_correction = atom->beta_correction;

	double force_multiplier = 0.;
	double friction_multiplier = 0.;
	double particle_energy = 0.;
	double exponent = 0.;
	double exp_cutoff = 20.;
	double one_plus_exp = 0.;
	double t_current;

    alpha_dot = 0.;

	double ndim = 3;

	if (nlocal != old_nlocal){
		memory->grow(f_corrected,nlocal,ndim,"fix:f_corrected");
		old_nlocal = nlocal;
	}

	//fprintf(screen,"pre assign final f_corrected [0,x] = %.9f \n",f_corrected[0][0]);
	//fprintf(screen,"pre assign final f [0,x] = %.9f \n",f);
	//fprintf(screen,"pre final beta_correction [0] = %.15f \n",beta_correction[0]);

	for(int i = 0; i < nlocal; ++i){
		for(int j = 0; j < ndim; ++j){
			f_corrected[i][j] = f[i][j];
		}
	}

	//fprintf(screen,"pre final calc f_corrected [0,x] = %.9f \n",f_corrected[0][0]);
	//fprintf(screen,"pre final calc f [0,x] = %.9f \n",f[0][0]);

	t_current = temperature->compute_scalar();
	//fprintf(screen,"t_current = %f \n",t_current);
	beta = 1 / (boltz * t_target);

	for(int i = 0; i < nlocal; ++i){
		if (mask[i] & groupbit){
			particle_energy = 0.;

			for(int j = 0; j < ndim; ++j){
				particle_energy += (v[i][j] * v[i][j]);
				//fprintf(screen,"final_integrate() v[i][j] = %d \n",v[i][j]);
			}
			particle_energy *= 0.5 * species_mass;
			//fprintf(screen,"initial_integrate() particle_energy = %f \n",particle_energy);
			exponent = beta * (particle_energy - mu);
			//fprintf(screen,"initial_integrate() exponent = %f \n",exponent);
			one_plus_exp = 1. + exp(exponent);
			//fprintf(screen,"initial_integrate() one_plus_exp = %f \n",one_plus_exp);

			if(exponent > exp_cutoff){
				alpha_dot += 0.5 * alpha_mass * species_mass * (2. * particle_energy * beta	-3.);

			} 
			else {
				alpha_dot -= 0.5 * alpha_mass * species_mass * (2. * particle_energy * beta * 
					(2. - one_plus_exp) / one_plus_exp + 3.) * (one_plus_exp - 1.) / one_plus_exp;
			}
		}
    }

	alpha_minus_quarter = alpha_plus_quarter;
	alpha_plus_quarter = alpha_minus_quarter + alpha_dot * dthalf;
	alpha_half = 0.5 * (alpha_minus_quarter + alpha_plus_quarter);

	//fprintf(screen,"final_integrate() alpha_half = %f \n",alpha_half);

  	for(int i = 0; i < nlocal; ++i){
  		if (mask[i] & groupbit){
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
				//apply beta_correction
				//fprintf(screen,"f_corrected[i][j] pre mod = %f \n",f_corrected[i][j]);
				f_corrected[i][j] = 1/(1-beta_correction[i]);
				//fprintf(screen,"f_corrected[i][j] post mod = %f \n",f_corrected[i][j]);

				v_minus_quarter[i][j] = v_plus_quarter[i][j];
		        potential_force = f_corrected[i][j] * force_multiplier;
		        friction_force = - friction_multiplier * alpha_half * v[i][j] * species_mass;
		        v_plus_quarter[i][j] = v_minus_quarter[i][j] + dthalf * (potential_force + friction_force) / species_mass;
				v[i][j] = 0.5 * (v_minus_quarter[i][j] + v_plus_quarter[i][j]);
	    	}
		}
  	}	

	//fprintf(screen,"post final calc f_corrected [0,x] = %.9f \n",f_corrected[0][0]);
	//fprintf(screen,"post final calc f [0,x] = %.9f \n",f[0][0]);
//  	t_current = temperature->compute_scalar();
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write - for fix_fermi the only variable
   that needs to be transmitted forward is alpha_plus_quarter.
------------------------------------------------------------------------- */

void FixFermi::write_restart(FILE *fp)
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

int FixFermi::pack_restart_data(double *list)
{
  int n = 0;
  list[n++] = alpha_plus_quarter;

  return n;
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixFermi::restart(char *buf)
{
  int n = 0;
  double *list = (double *) buf;
  alpha_plus_quarter = list[n++];
}
