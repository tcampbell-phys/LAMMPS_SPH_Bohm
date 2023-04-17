/* Bespoke Fix to assign electron velocities to a Fermi-Dirac kinetic energy distribution. Only
applicable to 3D systems with units electron.

Thomas Campbell */

#include "fix_init_fermi_vel.h"
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
#include "random_park.h"
#include <iostream>

using namespace LAMMPS_NS;
using namespace FixConst;

#define WARMUP 100

FixInitFermi::FixInitFermi(LAMMPS *lmp, int narg, char **arg):
	Fix(lmp, narg, arg), id_temp(NULL)
{
  if (narg < 4) error->all(FLERR,"Illegal fix InitFermi command");

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

  /* call system parameters for fermi distribution calculation*/

  t_target	= force->numeric(FLERR,arg[3]);
  //std::cout << "t_target =  " << t_target << std::endl;

  N_particle = force->numeric(FLERR,arg[4]);
  //std::cout << "N_particle =  " << N_particle << std::endl;
  
  mu = force->numeric(FLERR,arg[5]);
  //std::cout << "mu =  " << mu << std::endl;

  int seed = force->inumeric(FLERR,arg[6]);
  //std::cout << "seed =  " << seed << std::endl;

  double energy_lim = 5*mu;
  int N_int_points = 1*N_particle;
  int half_N_int_points = 0.5*N_int_points;
  double interval = energy_lim/N_int_points;

  //std::cout << "N_int_points =  " << N_int_points << std::endl;

  //std::cout << "interval =  " << interval << std::endl;

  memory->create(energy,N_int_points,"fix:energy");
  memory->create(cumul_dist,N_int_points,"fix:cumul_dist");
  memory->create(cumul_dist_norm,N_int_points,"fix:cumul_dist_norm");

  energy[0] = 0.0;
  cumul_dist[0] = 0.0;

  double** v = atom->v;

  int nlocal = atom->nlocal;
  int *mask = atom->mask;

  double rand_num;
  double distance;
  double raw_distance;
  double val_distance;
  double temp_distance;
  double particle_energy;
  double vx;
  double vy;
  double vz;
  double v_squared;
  double kin_E;
  double E_ratio;
  int jdx;

  //hard code in electron mass for now

  double e_mass = 0.0005485799;

  /* calculate cumulative Fermi distribution for given parameters */

  for(int i = 1; i < N_int_points; i++){
      energy[i] = energy[i-1] + interval;
      cumul_dist[i] = fermi_dist((energy[i]+energy[i-1])/2)*interval + cumul_dist[i-1];
  }

  for(int i = 1; i < N_int_points; i++){
      cumul_dist_norm[i] = cumul_dist[i]/cumul_dist[N_int_points-1];
  }

  /* assign particle energies using the cumulative Fermi distribution */

  int *type = atom->type;
  int dim = domain->dimension;

  int m;

  // create an atom map if one doesn't exist already

  int mapflag = 0;
  if (atom->map_style == 0) {
      mapflag = 1;
      atom->nghost = 0;
      atom->map_init();
      atom->map_set();
  }
  // error check

  if (atom->natoms > MAXSMALLINT)
      error->all(FLERR,"Too big a problem to use fix_init_fermi_vel create loop all");
  if (atom->tag_enable == 0)
      error->all(FLERR,
                "Cannot use fix_init_fermi_vel loop all unless atoms have IDs");
  if (atom->tag_consecutive() == 0)
      error->all(FLERR,
                "Atom IDs must be consecutive for fix_init_fermi_vel create loop all");

  // loop over all atoms in system
  // generate RNGs for all atoms, only assign to ones I own
  // use either per-type mass or per-atom rmass
  RanPark *random = NULL;
  random = new RanPark(lmp,seed);
  int natoms = static_cast<int> (atom->natoms);
  
  for(int i = 0; i <+ natoms; ++i){
    rand_num = random->uniform();
    raw_distance = cumul_dist_norm[0]-rand_num;
    distance = fabs(cumul_dist_norm[0]-rand_num);
    temp_distance = fabs(cumul_dist_norm[0]-rand_num);
    jdx = 0;
    for(int j = 1; j < N_int_points; j++){
      val_distance = fabs(cumul_dist_norm[j] - rand_num);
      if (val_distance < distance){
          temp_distance = fabs(cumul_dist_norm[j] - rand_num);
          distance = temp_distance;
          jdx = j;
      }
      if (val_distance > distance){
          particle_energy = energy[jdx];
          vx = random->uniform() - 0.5;
          vy = random->uniform() - 0.5;
          vz = random->uniform() - 0.5;

          v_squared = pow(vx,2) + pow(vy,2) + pow(vz,2);
          kin_E = 0.5 * e_mass * v_squared;

          E_ratio = kin_E/particle_energy;

          vx *= 1/(pow(E_ratio,0.5));
          vy *= 1/(pow(E_ratio,0.5));
          vz *= 1/(pow(E_ratio,0.5));

          m = atom->map(i);
          if (m >= 0 && m < nlocal) {
            if (mask[m] & groupbit) {
              v[m][0] = vx;
              v[m][1] = vy;
              v[m][2] = vz;
            }
          }
          break;
      }
    }
  }

  if (mapflag) {
    atom->map_delete();
    atom->map_style = 0;
    }
  delete random;
}

FixInitFermi::~FixInitFermi()
{
  memory->destroy(energy);
  memory->destroy(cumul_dist);
  memory->destroy(cumul_dist_norm);
}

int FixInitFermi::setmask()
{
  int mask = 0;
  return mask;
}

double FixInitFermi::fermi_dist(double E)
{
    boltz = 3.16681e-6;
    return (1/pow(M_PI,2))*pow(E,0.5)/(exp((E-mu)/(boltz*t_target))+1);
}
