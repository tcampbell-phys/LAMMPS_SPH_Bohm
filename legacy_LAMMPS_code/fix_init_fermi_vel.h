#ifdef FIX_CLASS
FixStyle(init_fermi_vel,FixInitFermi)
#else
#ifndef LMP_FIX_INITFERMI_H
#define LMP_FIX_INITFERMI_H

#include "fix.h"

namespace LAMMPS_NS {
class FixInitFermi : public Fix {
 public:
  FixInitFermi(class LAMMPS *, int, char **);
  virtual ~FixInitFermi();
  int setmask();

 protected:
  double fermi_dist(double E);
  double t_target;
  double mu;
  int N_particle;
  double boltz;
  double *energy;
  double *cumul_dist;
  double *cumul_dist_norm;

  char *id_temp;
};

}

#endif
#endif
