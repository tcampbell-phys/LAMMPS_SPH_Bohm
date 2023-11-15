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

PairStyle(coul/cut/SPH,PairCoulCutSPH)

#else

#ifndef LMP_PAIR_COUL_CUT_SPH_H
#define LMP_PAIR_COUL_CUT_SPH_H

#include "pair.h"

namespace LAMMPS_NS {

class PairCoulCutSPH : public Pair {
 public:
  PairCoulCutSPH(class LAMMPS *);
  virtual ~PairCoulCutSPH();
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
  int pseudo_key; // label to select pseudo parameters below
  int ion_charge; // ion charge
  double **cut,**scale;

  double sqrt2 = 1.4142135623730951;
  double sqrt_pi = 1.7724538509055159;

  //per-atom arrays
  double *theta_coul;
  double *chi_coul_ei;
  int *count;

  virtual void allocate();

  // Pseudopotential Gaussian Decomposition Parameters

// A - erf contribution removed
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

// B - no erf contribution removed. Slight deviation from true pseudo near origin
  int al_lda_B_key = 1;
  int al_lda_B_ion_charge = 3;
  int al_lda_B_Ncoeff = 7;
  double al_lda_B_a[7] = {0.002136611701287025, 0.11302714794414181, 0.023864321399827804, 16.358158351362377, 
  0.5303677735320211, 1.6243022290451936, 1.6122771477857079};
  double al_lda_B_c[7] = {-0.32036271967553065, -0.8633169812604677, -0.4424901366994751, 0.2138221574664385, 
  -2.399742395954945, -330.45034973789006, 335.90637417469407};
  
// C - no erf contribution removed. Slight deviation from true pseudo near origin
  int al_lda_C_key = 2;
  int al_lda_C_ion_charge = 3;
  int al_lda_C_Ncoeff = 11;
  double al_lda_C_a[11] = {0.017564281470838217, 0.0016115232238058735, 0.07792173292246932, 4.7392655944826245, 4.227393511766229, 0.3220501081344026, 
  4.339771532840177, 3.696809831601508, 4.858882850164823, 5.3228633595043355, 6.998884968039757};
  double al_lda_C_c[11] = {-0.3700172550561547, -0.27756311414555057, -0.6715971385710873, -961548.953125, -537981.71875, -1.4732938283123076, 
  807287.46875, 14116.133911132812, 720389.671875, -42448.140625, 189.86865234375};

// D - no erf contribution removed. Even temper generated.
  int al_lda_D_key = 3;
  int al_lda_D_ion_charge = 3;
  int al_lda_D_Ncoeff = 12;
  double al_lda_D_a[12] = {28.480358684358013, 10.0, 3.5111917342151315, 1.2328467394420664, 0.43287612810830595, 0.15199110829529341, 0.05336699231206311, 
  0.018738174228603847, 0.006579332246575677, 0.002310129700083161, 0.0008111308307896868, 0.0002848035868435803};
  double al_lda_D_c[12] = {0.04349978739628568, 0.10379401908721775, -0.8523820192203857, 5.8202089735423215, -1.9413628560432699, -0.5402857616718393, 
  -0.4928796991589479, -0.1168662776472047, -0.603258166462183, 1.6422574818134308, -3.336928315460682, 1.858105517923832};
  
  // E - GGA no erf contribution removed. Even temper generated.
  int al_gga_E_key = 4;
  int al_gga_E_ion_charge = 3;
  int al_gga_E_Ncoeff = 12;
  double al_gga_E_a[12] = {28.480358684358013, 10.0, 3.5111917342151315, 1.2328467394420664, 0.43287612810830595, 0.15199110829529341, 0.05336699231206311, 
  0.018738174228603847, 0.006579332246575677, 0.002310129700083161, 0.0008111308307896868, 0.0002848035868435803};
  double al_gga_E_c[12] = {-0.018256399111123756, 0.08399442658992484, -0.38052393405814655, 5.506835116888396, -1.7454933633562177, -0.6022056387446355, 
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

E: Pair style coul/cut requires atom attribute q

The atom style defined does not have these attributes.

*/
