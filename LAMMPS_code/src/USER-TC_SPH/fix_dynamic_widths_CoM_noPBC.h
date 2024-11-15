#ifdef FIX_CLASS
FixStyle(dynamic_widths_com_nopbc,FixDynamicWidthsCoMNoPBC)
#else
#ifndef LMP_FIX_DYNAMIC_WIDTHS_COM_NOPBC_H
#define LMP_FIX_DYNAMIC_WIDTHS_COM_NOPBC_H

#include "fix.h"

namespace LAMMPS_NS {
class FixDynamicWidthsCoMNoPBC : public Fix {
 public:
  FixDynamicWidthsCoMNoPBC(class LAMMPS *, int, char **);
  virtual ~FixDynamicWidthsCoMNoPBC();
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

  int nmax;

  double constant;
  double N_iter;
  double mix_fact;
  double start_width;
  double cut_global;
  double box_len;
  int type_avoid;
  int N_electrons;
  int N_elements_per_electron;
  int tag_ele_start;
  double cutsquared;
  double commflag;
  double assignment_flag;
  char *pair_name;

  double two_pi_over_len;

  double *x_mu;
  double *x_tau;
  double *y_mu;
  double *y_tau;
  double *z_mu;
  double *z_tau;

  double *mu,*tau;
  double *mu_all,*tau_all;

  int allocated;
  int unallocated;

  virtual void allocate();
  void deallocate();

  void FixedPointIterator();
  void CoM_Calculator();
  double Gauss_Width_Deriv(double pre_fact, double wid, double sep_sq);
 
};

}

#endif
#endif
