#ifdef FIX_CLASS
FixStyle(fixed_widths,FixFixedWidths)
#else
#ifndef LMP_FIX_FIXED_WIDTHS_H
#define LMP_FIX_FIXED_WIDTHS_H

#include "fix.h"

namespace LAMMPS_NS {
class FixFixedWidths : public Fix {
 public:
  FixFixedWidths(class LAMMPS *, int, char **);
  virtual ~FixFixedWidths();
  int setmask();
  virtual void init();
  virtual void setup_pre_force(int);
  virtual void setup_post_neighbor();
  virtual void min_pre_force(int);
  virtual void pre_force(int);
  virtual void post_neighbor();
  virtual int pack_forward_comm(int, int *, double *, int, int *);
  virtual void unpack_forward_comm(int, int, double *);
  virtual int pack_reverse_comm(int, int, double *);
  virtual void unpack_reverse_comm(int, int *, double *);

 protected:
  class Pair *pair;
  class NeighList *list;
  double constant;
  double commflag;
  double assignment_flag;
  double global_width;
  double cut_global;
  int type_avoid;
  double cutsquared;
  char *pair_name;
  void FixedPointIterator();
  double Gauss_Width_Deriv(double pre_fact, double wid, double sep_sq);
 
};

}

#endif
#endif
