/* Bespoke Fix to assign electron velocities to a Fermi-Dirac kinetic energy distribution. Only
applicable to 3D systems with units electron.

Thomas Campbell */

#include "fix_spherical_init.h"
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

FixSphericalInit::FixSphericalInit(LAMMPS *lmp, int narg, char **arg):
	Fix(lmp, narg, arg)
{
  if (narg < 5) error->all(FLERR,"Illegal fix SphericalInit command");

  /* call system parameters for fermi distribution calculation*/

  // Gaussian distribution h
  gauss_width	= force->numeric(FLERR,arg[3]);
  std::cout << "gauss_width =  " << gauss_width << std::endl;

  box_centre	= force->numeric(FLERR,arg[4]);
  std::cout << "box_centre =  " << box_centre << std::endl;

  N_particle = force->numeric(FLERR,arg[5]);

  int seed = force->inumeric(FLERR,arg[6]);
  std::cout << "seed =  " << seed << std::endl;

  double** x = atom->x;

  int nlocal = atom->nlocal;
  int *mask = atom->mask;

  double rand_num,rand_numb,rand_numc;
  double r_assign;
  double theta_assign;
  double phi_assign;

  double distance;
  double raw_distance;
  double val_distance;
  double temp_distance;
  int jdx;

  double radius_lim = 5*gauss_width;
  int N_int_points = 1*N_particle;
  int half_N_int_points = 0.5*N_int_points;
  double interval = radius_lim/N_int_points;

  memory->create(radius,N_int_points,"fix:radius");
  memory->create(cumul_dist,N_int_points,"fix:cumul_dist");
  memory->create(cumul_dist_norm,N_int_points,"fix:cumul_dist_norm");

  radius[0] = 0.0;
  cumul_dist[0] = 0.0;

  /* calculate cumulative Fermi distribution for given parameters */

  for(int i = 1; i < N_int_points; i++){
    radius[i] = radius[i-1] + interval;
    cumul_dist[i] = 2*gaussian_dens((radius[i]+radius[i-1])/2, gauss_width)*interval + cumul_dist[i-1];
  }

  fprintf(screen,"\n\ncumul_dist[-1] = %f",cumul_dist[N_int_points-1]);

  for(int i = 1; i < N_int_points; i++){
    cumul_dist_norm[i] = cumul_dist[i]/cumul_dist[N_int_points-1];
  }

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
      error->all(FLERR,"Too big a problem to use fix_spherical_init create loop all");
  if (atom->tag_enable == 0)
      error->all(FLERR,
                "Cannot use fix_spherical_init loop all unless atoms have IDs");
  if (atom->tag_consecutive() == 0)
      error->all(FLERR,
                "Atom IDs must be consecutive for fix_spherical_init create loop all");

  // loop over all atoms in system
  // generate RNGs for all atoms, only assign to ones I own
  // use either per-type mass or per-atom rmass
  RanPark *random = NULL;
  random = new RanPark(lmp,seed);
  int natoms = static_cast<int> (atom->natoms);
  
  for(int i = 0; i <+ natoms; ++i){
    m = atom->map(i);
    if (m >= 0 && m < nlocal) {
      if (mask[m] & groupbit) {

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
            r_assign = radius[jdx];
            break;
          }
        }

        rand_numb = random->uniform();
        theta_assign = M_PI*rand_numb;
        rand_numc = random->uniform();
        phi_assign = 2*M_PI*rand_numc;

        // fprintf(screen,"\n\nr_assign = %f",r_assign);
        // fprintf(screen,"\ntheta_assign = %f",theta_assign);
        // fprintf(screen,"\nphi_assign = %f",phi_assign);

        x[m][0] = box_centre + r_assign*sin(theta_assign)*cos(phi_assign);
        x[m][1] = box_centre + r_assign*sin(theta_assign)*sin(phi_assign);
        x[m][2] = box_centre + r_assign*cos(theta_assign);
        // fprintf(screen,"\nx[%d][0] = %f",m,x[m][0]);
        // fprintf(screen,"\nx[%d][1] = %f",m,x[m][1]);
        // fprintf(screen,"\nx[%d][2] = %f",m,x[m][2]);
      }
    }
  }

  if (mapflag) {
    atom->map_delete();
    atom->map_style = 0;
    }
  delete random;
}

FixSphericalInit::~FixSphericalInit()
{
  memory->destroy(radius);
  memory->destroy(cumul_dist);
  memory->destroy(cumul_dist_norm);
}

int FixSphericalInit::setmask()
{
  int mask = 0;
  return mask;
}

double FixSphericalInit::gaussian_dens(double radius_in, double width)
{
  return (1/(pow(2*M_PI,0.5)*width))*(exp(-(radius_in*radius_in)/(2*(width*width))));
}
