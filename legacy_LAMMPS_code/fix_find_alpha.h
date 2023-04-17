#ifdef FIX_CLASS
FixStyle(find_alpha,FixFindAlpha)
#else
#ifndef LMP_FIX_FINDALPHA_H
#define LMP_FIX_FINDALPHA_H

#include "fix.h"

namespace LAMMPS_NS {
class FixFindAlpha : public Fix {
 public:
  FixFindAlpha(class LAMMPS *, int, char **);
  int setmask();

 protected:
  double mu;
  double t_target;
  double N_particle;
  double boltz;
  char *id_temp;
};

}

#endif
#endif
