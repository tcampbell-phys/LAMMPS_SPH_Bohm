#ifdef FIX_CLASS
FixStyle(bohm_pot_calc_Mocz_rob_ftwo,FixBohmPotCalcMoczRobFtwo)
#else
#ifndef LMP_FIX_BOHM_POT_CALC_MOCZ_ROB_FTWO_H
#define LMP_FIX_BOHM_POT_CALC_MOCZ_ROB_FTWO_H

#include "fix.h"

namespace LAMMPS_NS {
class FixBohmPotCalcMoczRobFtwo : public Fix {
 public:
  FixBohmPotCalcMoczRobFtwo(class LAMMPS *, int, char **);
  virtual ~FixBohmPotCalcMoczRobFtwo();
  int setmask();
  virtual void init();
  virtual void end_of_step();
  virtual int pack_forward_comm(int, int *, double *, int, int *);
  virtual void unpack_forward_comm(int, int, double *);
  virtual int pack_reverse_comm(int, int, double *);
  virtual void unpack_reverse_comm(int, int *, double *);

 protected:
  class Pair *pair;
  class NeighList *list;
  int nmax;                         // allocated size of per-atom arrays

  // electron rest mass in AMU
  double e_mass = 0.00054857991;

  double cut_global;
  double gamma_factor;
  double hplanck;
  double hbar;
  double f_prefactor;
  char *pair_name;
  double cutsquared;
  double pi_fact;
  double **cut;

  int *numforce;

  // per-atom arrays

  double *dxx_rho_fix;
  double *dxy_rho_fix;
  double *dxz_rho_fix;
  double *dyy_rho_fix;
  double *dyz_rho_fix;
  double *dzz_rho_fix;

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

  int allocated;

  int commflag;
  void CalcBohmPot();
  
};

}

#endif
#endif
