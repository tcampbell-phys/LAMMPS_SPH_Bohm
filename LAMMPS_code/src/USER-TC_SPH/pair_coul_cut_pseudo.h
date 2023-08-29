/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(coul/cut/pseudo,PairCoulCutPseudo)

#else

#ifndef LMP_PAIR_COUL_CUT_PSEUDO_H
#define LMP_PAIR_COUL_CUT_PSEUDO_H

#include "pair.h"

namespace LAMMPS_NS {

class PairCoulCutPseudo : public Pair {
 public:
  PairCoulCutPseudo(class LAMMPS *);
  virtual ~PairCoulCutPseudo();
  virtual void compute(int, int);
  virtual void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  virtual void write_restart_settings(FILE *);
  virtual void read_restart_settings(FILE *);
  virtual double single(int, int, int, int, double, double, double, double &);
  void *extract(const char *, int &);

 protected:
  double cut_global;
  double **cut,**scale;
  int ion_species; // label which species are point charges
  int pseudo_key; // label to select pseudo parameters below
  int ion_charge; // ion charge

  virtual void allocate();

  // D - no erf contribution removed. Even temper generated.
  int al_lda_D_key = 3;
  int al_lda_D_ion_charge = 3;
  int al_lda_D_Ncoeff = 12;
  double al_lda_D_a[12] = {28.480358684358013, 10.0, 3.5111917342151315, 1.2328467394420664, 0.43287612810830595, 0.15199110829529341, 0.05336699231206311, 
  0.018738174228603847, 0.006579332246575677, 0.002310129700083161, 0.0008111308307896868, 0.0002848035868435803};
  double al_lda_D_c[12] = {0.04349978739628568, 0.10379401908721775, -0.8523820192203857, 5.8202089735423215, -1.9413628560432699, -0.5402857616718393, 
  -0.4928796991589479, -0.1168662776472047, -0.603258166462183, 1.6422574818134308, -3.336928315460682, 1.858105517923832};
  
  double *c_coeff;
  double *a_coeff;
  int N_coeff;
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
