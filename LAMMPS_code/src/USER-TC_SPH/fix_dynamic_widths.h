#ifdef FIX_CLASS
FixStyle(dynamic_widths,FixDynamicWidths)
#else
#ifndef LMP_FIX_DYNAMIC_WIDTHS_H
#define LMP_FIX_DYNAMIC_WIDTHS_H

#include "fix.h"
// #include "force.h"

namespace LAMMPS_NS {
class FixDynamicWidths : public Fix {
 public:
  FixDynamicWidths(class LAMMPS *, int, char **);
  virtual ~FixDynamicWidths();
  int setmask();
  void init();
  void post_integrate();

 protected:
  class Pair *pair;
  class NeighList *list; 
  // class PairHybrid *hybrid_pair;
  double constant;
  double N_iter;
  double mix_fact;
  double start_width;
  double cut_global;
  double cutsquared;
  char *pair_name;
  double Gauss_Width_Deriv(double pre_fact, double wid, double sep_sq);
 
};

}

#endif
#endif
