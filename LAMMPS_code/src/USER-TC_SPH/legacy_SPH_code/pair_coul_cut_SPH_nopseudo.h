/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Edited by Thomas Campbell (Oxford)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(coul/cut/SPH_nopseudo,PairCoulCutSPHNoPseudo)

#else

#ifndef LMP_PAIR_COUL_CUT_SPH_NOPSEUDO_H
#define LMP_PAIR_COUL_CUT_SPH_NOPSEUDO_H

#include "pair.h"

namespace LAMMPS_NS {

class PairCoulCutSPHNoPseudo : public Pair {
 public:
  PairCoulCutSPHNoPseudo(class LAMMPS *);
  virtual ~PairCoulCutSPHNoPseudo();
  virtual void compute(int, int);
  virtual void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
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
  int nmax;
  double cut_global;
  double ke_in;
  int ion_species; // label which species are point charges
  double **cut,**scale;

  double sqrt2 = 1.4142135623730951;
  double sqrt_pi = 1.7724538509055159;

  //per-atom arrays
  double *theta_coul;
  double *theta_coul_ei;
  int *count;

  virtual void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair style coul/cut requires atom attribute q

The atom style defined does not have these attributes.

*/
