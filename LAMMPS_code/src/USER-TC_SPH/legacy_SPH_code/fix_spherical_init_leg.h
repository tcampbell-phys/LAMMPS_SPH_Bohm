#ifdef FIX_CLASS
FixStyle(spherical_init,FixSphericalInit)
#else
#ifndef LMP_FIX_SPHERICAL_INIT_H
#define LMP_FIX_SPHERICAL_INIT_H

#include "fix.h"

namespace LAMMPS_NS {
class FixSphericalInit : public Fix {
 public:
  FixSphericalInit(class LAMMPS *, int, char **);
  virtual ~FixSphericalInit();
  int setmask();
  virtual void init();
  virtual void post_integrate();
  // virtual void setup_pre_force(int);
  // virtual void pre_force(int);
  // virtual void compute(int, int);
  virtual int pack_forward_comm(int, int *, double *, int, int *);
  virtual void unpack_forward_comm(int, int, double *);

 protected:
  class Pair *pair;
  class NeighList *list;
  double constant;
  double commflag;
  double cut_global;
  int type_avoid;
  int seed;
  double cutsquared;
  char *pair_name;
  void SphericalDistribution(double assign_radius_cutoff);
 
};

}

#endif
#endif
