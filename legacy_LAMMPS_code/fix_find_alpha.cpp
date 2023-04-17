/* Bespoke Fix to calculate alpha values (for electrons under the Fermi Nose-Hoover thermostat)
that return standard equations of motion.

Thomas Campbell */

#include "fix_find_alpha.h"
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
#include <mpi.h>

using namespace LAMMPS_NS;
using namespace FixConst;


FixFindAlpha::FixFindAlpha(LAMMPS *lmp, int narg, char **arg):
	Fix(lmp, narg, arg), id_temp(NULL)
{
  if (narg < 2) error->all(FLERR,"Illegal fix FindAlpha command");

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

  mu = force->numeric(FLERR,arg[3]);
  t_target = force->numeric(FLERR,arg[4]);
  N_particle = force->numeric(FLERR,arg[5]);

  double** v = atom->v;
  double** f = atom->f;

  int nlocal = atom->nlocal;
  int *mask = atom->mask;

  //hard code in electron mass for now

  double e_mass = 0.0005485799;

  double alpha = 0.0;
  double alpha_contrib;
  double v_squared;
  double force_vel_product;
  double zeta;
  double beta;
  double average_alpha;
  double E;
  double total_alpha;
  boltz = 3.16681e-6;
  beta = 1 / (boltz * t_target);

  for (int i = 0; i < nlocal; i++){
    if (mask[i] & groupbit){
      v_squared = pow(v[i][0],2) + pow(v[i][1],2) + pow(v[i][2],2);

      //std::cout << "i = "<< i << " v_squared = " << v_squared << std::endl;

      E = (0.5 * e_mass * v_squared) - mu;

      //std::cout << "i = "<< i << " E = " << E << std::endl;

      force_vel_product = f[i][0]*v[i][0] + f[i][1]*v[i][1] + f[i][2]*v[i][2];

      //std::cout << "i = "<< i << " force_vel_product = " << force_vel_product << std::endl;

      zeta = 1+exp(beta*E);

      //std::cout << "i = "<< i << " zeta = " << zeta << std::endl;

      alpha_contrib = ((2*force_vel_product*zeta)/(beta*pow(e_mass,2)*v_squared*exp(beta*E)))*(zeta*(log(zeta)-(beta*E))-1);

      alpha += alpha_contrib;

      std::cout << "i = "<< i << " alpha = " << alpha_contrib << std::endl;
    }
  }
  MPI_Allreduce(&alpha,&total_alpha,1,MPI_DOUBLE,MPI_SUM,world);
  average_alpha = total_alpha/N_particle; 
  std::cout << "ID: " << newarg[1] << " Calculated total alpha =" << total_alpha << std::endl;
  std::cout << "ID: " << newarg[1] << " Calculated average alpha =" << average_alpha << std::endl;
}

int FixFindAlpha::setmask()
{
  int mask = 0;
  return mask;
}

