/* Bespoke Fix to assign electron velocities to a Fermi-Dirac kinetic energy distribution. Only
applicable to 3D systems with units electron.

Thomas Campbell */

#include "fix_midrun_sample_vel.h"
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

FixMidrunSample::FixMidrunSample(LAMMPS *lmp, int narg, char **arg):
	Fix(lmp, narg, arg), id_temp(NULL)
{
  if (narg < 3) error->all(FLERR,"Illegal fix MidrunSample command");

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

  Nepe	= force->numeric(FLERR,arg[3]);
  std::cout << "Calling fix_midrun_sample_vel 1 with Nepe =   " << Nepe << std::endl;

  scaling_val	= force->numeric(FLERR,arg[4]);
  std::cout << "Calling fix_midrun_sample_vel 1 with scaling_val =   " << scaling_val << std::endl;

  int seed = force->inumeric(FLERR,arg[5]);

  double** v = atom->v;
  int *tagid = atom->tag;

  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  
  double mvv2e = force->mvv2e;

  double rand_num;
  double distance;
  double raw_distance;
  double val_distance;
  double temp_distance;
  double particle_energy;
  double vx_load;
  double vy_load;
  double vz_load;
  double vx;
  double vy;
  double vz;
  double v_squared;
  double kin_E;
  double E_ratio;
  int jdx;

  //hard code in electron mass for now

  double e_mass = 0.0005485799;

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
      error->all(FLERR,"Too big a problem to use fix_midrun_sample_vel create loop all");
  if (atom->tag_enable == 0)
      error->all(FLERR,
                "Cannot use fix_midrun_sample_vel loop all unless atoms have IDs");
  if (atom->tag_consecutive() == 0)
      error->all(FLERR,
                "Atom IDs must be consecutive for fix_midrun_sample_vel create loop all");

  // loop over all atoms in system
  // generate RNGs for all atoms, only assign to ones I own
  // use either per-type mass or per-atom rmass
  RanPark *random = NULL;
  random = new RanPark(lmp,seed);
  int natoms = static_cast<int> (atom->natoms);

  // // fprintf(screen,"\nnatoms = %d",natoms);
  
  for(int i = 0; i <+ natoms; ++i){
    rand_num = random->uniform();
    raw_distance = cumul_dist_norm[0]-rand_num;
    distance = fabs(cumul_dist_norm[0]-rand_num);
    temp_distance = fabs(cumul_dist_norm[0]-rand_num);
    jdx = 0;
    // fprintf(screen,"\ni = %d",i);
    
    for(int j = 1; j < N_int_points; j++){
      val_distance = fabs(cumul_dist_norm[j] - rand_num);
      // // fprintf(screen,"\nj = %d",j);
      if (val_distance < distance){
          temp_distance = fabs(cumul_dist_norm[j] - rand_num);
          distance = temp_distance;
          jdx = j;
      }
      if (val_distance > distance){
          // fprintf(screen,"\nfound val");
          vx = random->uniform() - 0.5;
          vy = random->uniform() - 0.5;
          vz = random->uniform() - 0.5;
          m = atom->map(i);
          // fprintf(screen,"\nCheck A");

          if (m >= 0 && m < nlocal) {
            if (mask[m] & groupbit) {

              vx_load = v[m][0];
              vy_load = v[m][1];
              vz_load = v[m][2];

              // fprintf(screen,"\nCheck B");

              particle_energy = energy[jdx]/(pow(Nepe,scaling_val));

              // fprintf(screen,"\nCheck C");

              v_squared = pow(vx_load,2) + pow(vy_load,2) + pow(vz_load,2);

              // fprintf(screen,"\nCheck D");
              kin_E = mvv2e * 0.5 * v_squared* e_mass/Nepe;

              // fprintf(screen,"\nCheck E");

              E_ratio = kin_E/particle_energy;

              // fprintf(screen,"\nCheck F");
  
              // fprintf(screen,"\nCheck G");
              v[m][0] *= 1/(pow(E_ratio,0.5));
              v[m][1] *= 1/(pow(E_ratio,0.5));
              v[m][2] *= 1/(pow(E_ratio,0.5));
              // fprintf(screen,"\nCheck H");
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

FixMidrunSample::~FixMidrunSample()
{
  return;
}

int FixMidrunSample::setmask()
{
  int mask = 0;
  return mask;
}