/* -*- c++ -*- ----------------------------------------------------------
Thomas Campbell (Oxford)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(filinovs,PairFilinovSinglet)

#else

#ifndef LMP_PAIR_FILINOV_SINGLET_H
#define LMP_PAIR_FILINOV_SINGLET_H

#include "pair.h"
#include <vector>

namespace LAMMPS_NS {

class PairFilinovSinglet : public Pair {
 public:
  PairFilinovSinglet(class LAMMPS *);
  virtual ~PairFilinovSinglet();
  virtual void compute(int, int);
  virtual void settings(int, char **);
  void coeff(int, char **);
  virtual void init_style();
  virtual double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  virtual void write_restart_settings(FILE *);
  virtual void read_restart_settings(FILE *);
  // virtual double single(int, int, int, int, double, double, double, double &);
  virtual void *extract(const char *, int &);

 protected:
  double cut_global;
  double targ_temp; //Kelvin for units electron
  double boltz_val;
  double hbar_val;
  double e_mass;
  double **cut;

  virtual void allocate();
};

}

#endif
#endif
