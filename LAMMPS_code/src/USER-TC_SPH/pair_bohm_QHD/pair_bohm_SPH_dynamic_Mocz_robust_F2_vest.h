/* -*- c++ -*- ----------------------------------------------------------
Bespoke pair_style to compute the Bohm force using an SPH-style pressure
tensor, with dynamic per-particle widths.

pair_bohm_SPH_dynamic:
- Electrons as basic fluid element.
- Dynamic per-particle gaussian widths.
- Mocz 2015. form for second derivatives.
- Bohm potential calculated on total density of electron fluid (QHD Bohm Potential).

Thomas Campbell
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(bohm_SPH_dynamic_Mocz_rob_ftwo_vest,PairBohmSPHDynamicMoczRobustFtwoVest)

#else

#ifndef LMP_PAIR_BOHM_SPH_DYNAMIC_MOCZ_ROBUST_FTWO_VEST_H
#define LMP_PAIR_BOHM_SPH_DYNAMIC_MOCZ_ROBUST_FTWO_VEST_H

#include "pair.h"

namespace LAMMPS_NS {

class PairBohmSPHDynamicMoczRobustFtwoVest : public Pair {
 public:
  PairBohmSPHDynamicMoczRobustFtwoVest(class LAMMPS *);
  virtual ~PairBohmSPHDynamicMoczRobustFtwoVest();

  virtual void compute(int, int);
  virtual void settings(int, char **);
  void coeff(int, char **);
  virtual void init_style();
  virtual double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  virtual void write_restart_settings(FILE *);
  virtual void read_restart_settings(FILE *);
  virtual int pack_forward_comm(int, int *, double *, int, int *);
  virtual void unpack_forward_comm(int, int, double *);
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  void *extract(const char *, int &);
 

protected:
  int nmax;                         // allocated size of per-atom arrays

  // electron rest mass in AMU
  double e_mass = 0.00054857991;

  double cut_global;
  double gamma_factor;
  double hplanck;
  double hbar;
  double f_prefactor;
  double cutsquared;
  double pi_fact;
  double **cut;

  int *numforce;

  // per-atom arrays

  double *dxx_rho;
  double *dxy_rho;
  double *dxz_rho;
  double *dyy_rho;
  double *dyz_rho;
  double *dzz_rho;

  // compute dx_rho terms internally in Bohm algorithm as F2 derivatives

  double *dx_rho;
  double *dy_rho;
  double *dz_rho;

  double *Pxx;
  double *Pxy;
  double *Pxz;
  double *Pyz;
  double *Pyy;
  double *Pzz;


  int commflag;


  virtual void allocate();
};

}

#endif
#endif