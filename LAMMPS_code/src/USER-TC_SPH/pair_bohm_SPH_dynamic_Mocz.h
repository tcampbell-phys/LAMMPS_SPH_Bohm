/* -*- c++ -*- ----------------------------------------------------------
Bespoke pair_style to compute the Bohm force using an SPH-style pressure
tensor, with dynamic per-particle widths.

pair_bohm_SPH_dynamic:
- Electrons as basic fluid element.
- Dynamic per-particle gaussian widths.
- Mocz 2015. form for second derivatives.

Thomas Campbell
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(bohm_SPH_dynamic_Mocz,PairBohmSPHDynamicMocz)

#else

#ifndef LMP_PAIR_BOHM_SPH_DYNAMIC_MOCZ_H
#define LMP_PAIR_BOHM_SPH_DYNAMIC_MOCZ_H

#include "pair.h"

namespace LAMMPS_NS {

class PairBohmSPHDynamicMocz : public Pair {
 public:
  PairBohmSPHDynamicMocz(class LAMMPS *);
  virtual ~PairBohmSPHDynamicMocz();

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
  double e_mass = 0.0005485799;

  double cut_global;
  double gamma_factor;
  double hplanck;
  double hbar;
  double f_prefactor;
  double **cut;

  int *numforce;

  // per-atom arrays

  double *dxx_rho;
  double *dxy_rho;
  double *dxz_rho;
  double *dyy_rho;
  double *dyz_rho;
  double *dzz_rho;


  int commflag;


  virtual void allocate();
};

}

#endif
#endif