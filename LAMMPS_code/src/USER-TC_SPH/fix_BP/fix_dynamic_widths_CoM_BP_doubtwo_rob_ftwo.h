#ifdef FIX_CLASS
FixStyle(dynamic_widths_com_bp_doubtwo_rob_ftwo,FixDynamicWidthsCoMBPDoubtwoRobFtwo)
#else
#ifndef LMP_FIX_DYNAMIC_WIDTHS_COM_BP_DOUBTWO_ROB_FTWO_H
#define LMP_FIX_DYNAMIC_WIDTHS_COM_BP_DOUBTWO_ROB_FTWO_H

#include "fix.h"

namespace LAMMPS_NS {
class FixDynamicWidthsCoMBPDoubtwoRobFtwo : public Fix {
 public:
  FixDynamicWidthsCoMBPDoubtwoRobFtwo(class LAMMPS *, int, char **);
  virtual ~FixDynamicWidthsCoMBPDoubtwoRobFtwo();
  int setmask();
  virtual void init();
  virtual void end_of_step();
  virtual void setup_pre_force(int);
  virtual void setup_post_neighbor();
  virtual void min_pre_force(int);
  virtual void pre_force(int);
  virtual void post_neighbor();
  virtual int pack_forward_comm(int, int *, double *, int, int *);
  virtual void unpack_forward_comm(int, int, double *);
  virtual int pack_reverse_comm(int, int, double *);
  virtual void unpack_reverse_comm(int, int *, double *);
  double compute_scalar();

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
  double pi_fact;
  char *pair_name;

  // electron rest mass in AMU
  double e_mass = 0.00054857991;

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

  // Bohm Potential variables

  double gamma_factor;
  double hplanck;
  double hbar;
  double f_prefactor;

  // per-atom arrays

  double *dxx_rho_fix;
  double *dxy_rho_fix;
  double *dxz_rho_fix;
  double *dyy_rho_fix;
  double *dyz_rho_fix;
  double *dzz_rho_fix;
  double *dyx_rho_fix;
  double *dzx_rho_fix;
  double *dzy_rho_fix;

  // compute dx_rho_fix terms internally in Bohm algorithm as F2 derivatives

  double *dx_rho_fix;
  double *dy_rho_fix;
  double *dz_rho_fix;

  double *Pxx_fix;
  double *Pxy_fix;
  double *Pxz_fix;
  double *Pyz_fix;
  double *Pyy_fix;
  double *Pzz_fix;
  double *Pyx_fix;
  double *Pzx_fix;
  double *Pzy_fix;

  int eflag;
  double eBohm,eBohm_all;

  virtual void allocate();
  void deallocate();

  void FixedPointIterator();
  void CoM_Calculator();
  void CalcBohmPot();
  double Gauss_Width_Deriv(double pre_fact, double wid, double sep_sq);
 
};

}

#endif
#endif
