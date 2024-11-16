/* -*- c++ -*- ----------------------------------------------------------
Thomas Campbell (Oxford)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(confine_full,PairConfineFull)

#else

#ifndef LMP_PAIR_CONFINE_FULL_H
#define LMP_PAIR_CONFINE_FULL_H

#include "pair.h"
#include <vector>

namespace LAMMPS_NS {

class PairConfineFull : public Pair {
 public:
  PairConfineFull(class LAMMPS *);
  virtual ~PairConfineFull();
  virtual void compute(int, int);
  virtual void settings(int, char **);
  void coeff(int, char **);
  virtual void init_style();
  virtual double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  virtual void write_restart_settings(FILE *);
  virtual void read_restart_settings(FILE *);
  virtual void *extract(const char *, int &);

 protected:
  double strength;
  double box_len;
  double half_box_len;
  double N_epe;
  int N_epe_minone;
  int N_SPH;
  int tag_ele_start;

  double force_factor,force_add_factor;
  double cut_global;
  double boltz_val;
  double hbar_val;
  double **cut;

  double *all_delx;
  double *all_dely;
  double *all_delz;

  double *ind_delx;
  double *ind_dely;
  double *ind_delz;

  int unallocated;

  virtual void allocate();
  void deallocate();
};

}

#endif
#endif
