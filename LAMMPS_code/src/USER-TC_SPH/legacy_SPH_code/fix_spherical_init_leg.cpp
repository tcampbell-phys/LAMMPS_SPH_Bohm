/* -*- c++ -*- ----------------------------------------------------------
Bespoke fix to compute the updated with valeus via a fixed point iteration,
and associated omega_SPH values for use in SPH force expressions.

Thomas Campbell
------------------------------------------------------------------------- */

#include "fix_spherical_init.h"
#include <cstring>
#include <mpi.h>
#include "group.h"
#include "modify.h"
#include "error.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "pair.h"
#include "pair_hybrid.h"
#include <cmath>
#include <math.h>
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "irregular.h"
#include "compute.h"
#include "domain.h"
#include "utils.h"
#include "memory.h"
#include "random_park.h"

using namespace LAMMPS_NS;
using namespace FixConst;

FixSphericalInit::FixSphericalInit(LAMMPS *lmp, int narg, char **arg):
	Fix(lmp, narg, arg)
{
  if (narg < 6) error->all(FLERR,"Illegal fix Fixed Widths command");

  cut_global = force->numeric(FLERR,arg[3]);
  fprintf(screen,"\nCutoff spherical radius = %f",cut_global);
  // type not involved in SPH loop
  type_avoid = force->numeric(FLERR,arg[4]);
  seed = force->inumeric(FLERR,arg[5]);
  pair_name = strdup(arg[6]);

  cutsquared = cut_global*cut_global;

  comm_forward = 3;
  comm_reverse = 3; 

  commflag = 0;
}

FixSphericalInit::~FixSphericalInit()
{
  //free pair_name variable
  free(pair_name);
}

int FixSphericalInit::setmask()
{
  int mask = 0;
  // mask |= FixConst::PRE_FORCE;
  mask |= FixConst::POST_INTEGRATE;
  return mask;
}

void FixSphericalInit::init()
{ 
  double *width_SPH = atom->width_SPH;
  double *rho_SPH = atom->rho_SPH;
  double *omega_SPH = atom->omega_SPH;
  double *u_SPH = atom->u_SPH;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  pair = force->pair;
  // If a hybrid style is used we need to acces the correct sub-style.
  PairHybrid *hybrid_pair = dynamic_cast<PairHybrid*> (pair);

  if (hybrid_pair) {
    // The pair style is a hybrid style.
    if (!pair_name) error->all(FLERR,"When a hybrid pair-style is used, 'pair_name' must be set for the lagrangian solver.");
    int nstyles = hybrid_pair->nstyles;
    int found = -1;
    for (int m = 0; m < nstyles; ++m) {
      if (strcmp(pair_name, hybrid_pair->keywords[m]) == 0) {
        if (found != -1) error->all(FLERR,"Multiple pair-styles with the name 'pair_name' found.");
        found = m;
      }
    }
    if (found == -1) error->all(FLERR,"No pair-style with the name of 'pair_name' found.");
    // Set the correct pair style.
    pair = hybrid_pair->styles[found];
  }
}


// void FixSphericalInit::setup_pre_force(int)
// {
//   SphericalDistribution(cut_global);
// }

// void FixSphericalInit::pre_force(int)
// {
//  int ntimestep = update->ntimestep; // current timestep
//   if (ntimestep == 0) {
//     SphericalDistribution(cut_global);
//   }
// }

void FixSphericalInit::post_integrate()
{
 int ntimestep = update->ntimestep; // current timestep
  if (ntimestep == 1) {
    SphericalDistribution(cut_global);
  }
}

// void PairCoulCutSPHNoPseudoSubFix::compute(int eflag, int vflag)
// {
//   int ntimestep = update->ntimestep; // current timestep
//   if (ntimestep == 0) {
//     SphericalDistribution(cut_global);
//   }
// }

void FixSphericalInit::SphericalDistribution(double assign_radius_cutoff)
{
  int a,i,j,ii,jj,inum,jnum,itype,jtype;
  int *ilist,*jlist,*numneigh,**firstneigh;
	double **x = atom->x;

	int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int newton_pair = force->newton_pair;

  // spherical assignment procedure
  double rand_num,rand_numb,rand_numc;
  double r_assign;
  double theta_assign;
  double phi_assign;

  list = pair->list;
	
  inum = list->inum;
  ilist = list->ilist;

  RanPark *random = NULL;
  random = new RanPark(lmp,seed);

  fprintf(screen,"\n FixSphericalInit::SphericalDistribution() called...");

  for(ii = 0; ii < inum; ii++){
    i = ilist[ii];

    itype = type[i];
    if (type[i] == type_avoid){
      continue;
    }

    rand_num = random->uniform();
    r_assign = rand_num*assign_radius_cutoff;
    rand_numb = random->uniform();
    theta_assign = M_PI*rand_numb;
    rand_numc = random->uniform();
    phi_assign = 2*M_PI*rand_numc;

    fprintf(screen,"\n\nr_assign = %f",r_assign);
    fprintf(screen,"\ntheta_assign = %f",theta_assign);
    fprintf(screen,"\nphi_assign = %f",phi_assign);

    x[i][0] = r_assign*sin(theta_assign)*cos(phi_assign);
    x[i][1] = r_assign*sin(theta_assign)*sin(phi_assign);
    x[i][2] = r_assign*cos(theta_assign);
  }

  comm_forward = 3;
  comm_reverse = 3;

  comm->forward_comm_fix(this);
}


int FixSphericalInit::pack_forward_comm(int n, int *list, double *buf,
                               int /*pbc_flag*/, int * /*pbc*/)
{
  int i,j,m;

  m = 0;
  double **x = atom->x;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = x[i][0];
    buf[m++] = x[i][1];
    buf[m++] = x[i][2];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixSphericalInit::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  double **x = atom->x;
  for (i = first; i < last; i++){
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
  }
  
}