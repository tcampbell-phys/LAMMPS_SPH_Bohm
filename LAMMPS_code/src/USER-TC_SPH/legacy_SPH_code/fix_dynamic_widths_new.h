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
  virtual void init();
  virtual void setup_post_neighbor();
  virtual void pre_force(int);
  virtual void post_neighbor();
  virtual void min_pre_force(int);
  virtual int pack_forward_comm(int, int *, double *, int, int *);
  virtual void unpack_forward_comm(int, int, double *);
  virtual int pack_reverse_comm(int, int, double *);
  virtual void unpack_reverse_comm(int, int *, double *);

 protected:
  class Pair *pair;
  class NeighList *list;
  double constant;
  double N_iter;
  double mix_fact;
  double start_width;
  double cut_global;
  int type_avoid;
  double cutsquared;
  double commflag;
  char *pair_name;
  void FixedPointIterator();
  double Gauss_Width_Deriv(double pre_fact, double wid, double sep_sq);
 
};

}

#endif
#endif
