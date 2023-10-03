/* -*- c++ -*- ----------------------------------------------------------
Thomas Campbell (Oxford)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(filinovt,PairFilinovTriplet)

#else

#ifndef LMP_PAIR_FILINOV_TRIPLET_H
#define LMP_PAIR_FILINOV_TRIPLET_H

#include "pair.h"
#include <vector>

namespace LAMMPS_NS {

class PairFilinovTriplet : public Pair {
 public:
  PairFilinovTriplet(class LAMMPS *);
  virtual ~PairFilinovTriplet();
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
  double hh_me;
  double e_mass;
  double **cut;

  virtual void allocate();
};

}

#endif
#endif
