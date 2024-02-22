/* -*- c++ -*- ----------------------------------------------------------
Thomas Campbell (Oxford)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(confine_fixed,PairConfineFixed)

#else

#ifndef LMP_PAIR_CONFINE_FIXED_H
#define LMP_PAIR_CONFINE_FIXED_H

#include "pair.h"
#include <vector>

namespace LAMMPS_NS {

class PairConfineFixed : public Pair {
 public:
  PairConfineFixed(class LAMMPS *);
  virtual ~PairConfineFixed();
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

  double cut_global;
  double x_fix, y_fix, z_fix;
  double boltz_val;
  double hbar_val;
  double **cut;

  virtual void allocate();
};

}

#endif
#endif
