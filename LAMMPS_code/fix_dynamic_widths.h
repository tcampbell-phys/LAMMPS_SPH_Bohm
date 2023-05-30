#ifdef FIX_CLASS
FixStyle(dynamic_widths,FixDynamicWidths)
#else
#ifndef LMP_FIX_DYNAMIC_WIDTHS_H
#define LMP_FIX_DYNAMIC_WIDTHS_H

#include "fix.h"

namespace LAMMPS_NS {
class FixDynamicWidths : public Fix {
 public:
  FixDynamicWidths(class LAMMPS *, int, char **);
  virtual ~FixDynamicWidths();
  int setmask();
  void init();
  void min_pre_force();
  void pre_force();

 protected:
  double constant;
  double N_iter;
  double mix_fact;
  double start_width;
  double cut_global;
  double Gauss_Width_deriv(double wid, double sep);
 
}

}

#endif
#endif
