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

PairStyle(sub_coul/long/SPH_ph,PairSubCoulLongSPHph)

#else

#ifndef LMP_PAIR_COUL_LONG_SPH_PH_H
#define LMP_PAIR_COUL_LONG_SPH_PH_H

#include "pair.h"
#include <vector>

namespace LAMMPS_NS {

class PairSubCoulLongSPHph : public Pair {
 public:
  PairSubCoulLongSPHph(class LAMMPS *);
  ~PairSubCoulLongSPHph();
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
//   virtual double single(int, int, int, int, double, double, double, double &);
  virtual void *extract(const char *, int &);

 protected:
  int nmax; // allocated size of per-atom arrays

  double cut_coul,cut_coulsq,qdist;
  double *cut_respa;
  double g_ewald;
  double ke_in; //currently redundant - could be used as a multiplier
  int ion_species; // label which species are point charges
  int pseudo_key; // label to select pseudo parameters below
  int ion_charge; // ion charge
  double ele_TFWHM; // twice the full width half maximum of the average electron in SPH system
  double input_g_ewald;
  double **scale;

  // per-atom arrays
  double *theta_coul;
  double *chi_coul_ei;

  int numforce; // value to check whether many body quantities have been computed prior to calling single()

  double sqrt2 = 1.4142135623730951;
  double sqrt_pi = 1.7724538509055159;

  virtual void allocate();

  // Pseudopotential Gaussian Decomposition Parameters

  // A - al_lda for 5.2 g/cc 3+ ions, erf contribution removed
  int al_lda_A_key = 0;
  int al_lda_A_ion_charge = 3;
  int al_lda_A_Ncoeff = 13;
  double al_lda_A_a[13] = {0.003876418186049235, 1.150625641772063, 0.056097205005079455, 0.019689359264809383, 
  1.747983820726069, 1.0998362686176768, 1.305314808231124, 0.36011457237236427, 
  1.5545126182845386, 1.9257077692314741, 1.8243933672383936, 7.4009573484454805, 
  0.10506523220666016};
  double al_lda_A_c[13] = {-0.20478337125996404, -22542.599365234375, -0.29884432538528927, -0.3111565024664742, 
  89734.388671875, 11642.630920410156, 22648.13720703125, -2.8933730144053698, 
  -36925.607421875, 18654.174072265625, -83205.1767578125, -0.14233727008104324, 
  -0.41864406527020037};

  // D - LDA no erf contribution removed. Even temper generated.
  int al_lda_D_key = 3;
  int al_lda_D_ion_charge = 3;
  int al_lda_D_Ncoeff = 12;
  double al_lda_D_a[12] = {28.480358684358013, 10.0, 3.5111917342151315, 1.2328467394420664, 0.43287612810830595, 0.15199110829529341, 0.05336699231206311, 
  0.018738174228603847, 0.006579332246575677, 0.002310129700083161, 0.0008111308307896868, 0.0002848035868435803};
  double al_lda_D_c[12] = {0.04349978739628568, 0.10379401908721775, -0.8523820192203857, 5.8202089735423215, -1.9413628560432699, -0.5402857616718393, 
  -0.4928796991589479, -0.1168662776472047, -0.603258166462183, 1.6422574818134308, -3.336928315460682, 1.858105517923832};

  // E - GGA no erf contribution removed. Even temper generated.
  int al_gga_A_key = 4;
  int al_gga_A_ion_charge = 3;
  int al_gga_A_Ncoeff = 12;
  double al_gga_A_a[12] = {28.480358684358013, 10.0, 3.5111917342151315, 1.2328467394420664, 0.43287612810830595, 0.15199110829529341, 0.05336699231206311, 
  0.018738174228603847, 0.006579332246575677, 0.002310129700083161, 0.0008111308307896868, 0.0002848035868435803};
  double al_gga_A_c[12] = {-0.018256399111123756, 0.08399442658992484, -0.38052393405814655, 5.506835116888396, -1.7454933633562177, -0.6022056387446355, 
  -0.4560560590471141, -0.17463317525107414, -0.385358068626374, 0.8103848937898874, -1.733453780412674, 0.9091406837105751};
  
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

E: Pair style lj/cut/coul/long requires atom attribute q

The atom style defined does not have this attribute.

E: Pair style requires a KSpace style

No kspace style is defined.

*/
