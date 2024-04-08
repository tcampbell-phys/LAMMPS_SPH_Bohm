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

 protected:

  double gaussian_dens(double radius_in, double width);

  double gauss_width;
  double box_centre;
  double N_particle;

  double *radius;
  double *cumul_dist;
  double *cumul_dist_norm;
};

}

#endif
#endif
