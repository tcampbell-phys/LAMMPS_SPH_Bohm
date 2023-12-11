/* -*- c++ -*- ----------------------------------------------------------
Thomas Campbell (Oxford)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(confine,PairConfine)

#else

#ifndef LMP_PAIR_CONFINE_H
#define LMP_PAIR_CONFINE_H

#include "pair.h"
#include <vector>

namespace LAMMPS_NS {

class PairConfine : public Pair {
 public:
  PairConfine(class LAMMPS *);
  virtual ~PairConfine();
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
  double cut_global;
  double strength;
  double N_elements_per_electron;
  double **cut;

  int tag_ele_start;

  virtual void allocate();
};

}

#endif
#endif
