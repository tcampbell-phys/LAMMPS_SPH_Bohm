/* Fix to change Force computation in the Bohm GWP scheme to ensure energy conservation*/

#ifdef FIX_CLASS

FixStyle(conserve,FixConserve)

#else

#ifndef LMP_FIX_CONSERVE_H
#define LMP_FIX_CONSERVE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixConserve : public Fix {
 public:
  FixConserve(class LAMMPS *, int, char **);
  int setmask();
  void post_force();
};

}

#endif
#endif