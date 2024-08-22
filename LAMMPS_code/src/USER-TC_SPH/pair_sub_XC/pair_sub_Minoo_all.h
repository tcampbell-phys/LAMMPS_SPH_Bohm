/* -*- c++ -*- ----------------------------------------------------------
Thomas Campbell (Oxford)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(sub_minoo_all,PairSubMinooAll)

#else

#ifndef LMP_PAIR_SUB_MINOO_ALL_H
#define LMP_PAIR_SUB_MINOO_ALL_H

#include "pair.h"
#include <vector>

namespace LAMMPS_NS {

class PairSubMinooAll : public Pair {
 public:
  PairSubMinooAll(class LAMMPS *);
  virtual ~PairSubMinooAll();
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
  double N_elements_per_electron;
  double N_elements_per_electron_sq;
  int tag_ele_start;
  double hbar_val;
  double hh_me;
  double e_mass;
  double pot_pre_fact;
  double for_pre_fact;
  double exp_fact;

  double ln_two = 0.6931471805599453;
  double **cut;

  virtual void allocate();
};

}

#endif
#endif
