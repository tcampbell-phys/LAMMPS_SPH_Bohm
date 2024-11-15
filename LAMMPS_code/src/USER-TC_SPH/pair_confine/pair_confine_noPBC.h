/* -*- c++ -*- ----------------------------------------------------------
Thomas Campbell (Oxford)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(confine_nopbc,PairConfineNoPBC)

#else

#ifndef LMP_PAIR_CONFINE_NOPBC_H
#define LMP_PAIR_CONFINE_NOPBC_H

#include "pair.h"
#include <vector>

namespace LAMMPS_NS {

class PairConfineNoPBC : public Pair {
 public:
  PairConfineNoPBC(class LAMMPS *);
  virtual ~PairConfineNoPBC();
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

  double force_factor;
  double cut_global;
  double boltz_val;
  double hbar_val;
  double **cut;

  virtual void allocate();
};

}

#endif
#endif
