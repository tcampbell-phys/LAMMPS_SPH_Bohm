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

PairStyle(coul/long/SPH_ph,PairCoulLongSPHph)

#else

#ifndef LMP_PAIR_COUL_LONG_SPH_PH_H
#define LMP_PAIR_COUL_LONG_SPH_PH_H

#include "pair.h"
#include <vector>

namespace LAMMPS_NS {

class PairCoulLongSPHph : public Pair {
 public:
  PairCoulLongSPHph(class LAMMPS *);
  ~PairCoulLongSPHph();
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
  int N_elements_per_electron; //SPH elements per electron
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
  int al_gga_E_key = 4;
  int al_gga_E_ion_charge = 3;
  int al_gga_E_Ncoeff = 12;
  double al_gga_E_a[12] = {28.480358684358013, 10.0, 3.5111917342151315, 1.2328467394420664, 0.43287612810830595, 0.15199110829529341, 0.05336699231206311, 
  0.018738174228603847, 0.006579332246575677, 0.002310129700083161, 0.0008111308307896868, 0.0002848035868435803};
  double al_gga_E_c[12] = {-0.018256399111123756, 0.08399442658992484, -0.38052393405814655, 5.506835116888396, -1.7454933633562177, -0.6022056387446355, 
  -0.4560560590471141, -0.17463317525107414, -0.385358068626374, 0.8103848937898874, -1.733453780412674, 0.9091406837105751};

  // F - scaled (for 10 elements per electron) GGA no erf contribution removed. Even temper generated.
  int al_gga_F_key = 100;
  int al_gga_F_ion_charge = 3;
  int al_gga_F_Ncoeff = 12;
  int al_gga_F_Nepe = 10;
  double al_gga_F_a[12] = {28.480358684358013, 10.0, 3.5111917342151315, 1.2328467394420664, 0.43287612810830595, 0.15199110829529341, 0.05336699231206311, 
  0.018738174228603847, 0.006579332246575677, 0.002310129700083161, 0.0008111308307896868, 0.0002848035868435803};
  double al_gga_F_c[12] = {-0.0018256399162055459, -0.03805239340726985, 0.008399442667723633, 0.5506835116975708, -0.17454933634144254, -0.06022056386063923, -0.04560560591562535, 
  -0.017463317504734732, -0.038535806990694255,  0.08103849040344357, -0.17334537860006094, 0.09091406920924783};

  // G - scaled (for 8 elements per electron) GGA no erf contribution removed. Even temper generated.
  int al_gga_G_key = 101;
  int al_gga_G_ion_charge = 3;
  int al_gga_G_Ncoeff = 12;
  int al_gga_G_Nepe = 8;
  double al_gga_G_a[12] = {28.480358684358013, 10.0, 3.5111917342151315, 1.2328467394420664, 0.43287612810830595, 0.15199110829529341, 0.05336699231206311, 
  0.018738174228603847, 0.006579332246575677, 0.002310129700083161, 0.0008111308307896868, 0.0002848035868435803};
  double al_gga_G_c[12] = {-0.0022820498888904694, 0.010499303323740605, -0.04756549175726832, 0.6883543896110496, -0.21818667041952722, -0.07527570484307944, 
  -0.05700700738088926, -0.021829146906384267, -0.04816975857829675, 0.10129811172373593, -0.21668172255158424, 0.11364258546382189};
  
  // H - scaled (for 4 elements per electron) GGA no erf contribution removed. Even temper generated.
  int al_gga_H_key = 102;
  int al_gga_H_ion_charge = 3;
  int al_gga_H_Ncoeff = 12;
  int al_gga_H_Nepe = 4;
  double al_gga_H_a[12] = {28.480358684358013, 10.0, 3.5111917342151315, 1.2328467394420664, 0.43287612810830595, 0.15199110829529341, 0.05336699231206311, 
  0.018738174228603847, 0.006579332246575677, 0.002310129700083161, 0.0008111308307896868, 0.0002848035868435803};
  double al_gga_H_c[12] = {-0.004564099777780939, 0.02099860664748121, -0.09513098351453664, 1.376708779222099, -0.43637334083905444, -0.15055140968615888, -0.11401401476177853, 
  -0.043658293812768534, -0.0963395171565935, 0.20259622344747186, -0.4333634451031685, 0.22728517092764378};

  // I - scaled (for 15 elements per electron) GGA no erf contribution removed. Even temper generated.
  int al_gga_I_key = 103;
  int al_gga_I_ion_charge = 3;
  int al_gga_I_Ncoeff = 12;
  int al_gga_I_Nepe = 15;
  double al_gga_I_a[12] = {28.480358684358013, 10.0, 3.5111917342151315, 1.2328467394420664, 0.43287612810830595, 0.15199110829529341, 0.05336699231206311, 
  0.018738174228603847, 0.006579332246575677, 0.002310129700083161, 0.0008111308307896868, 0.0002848035868435803};
  double al_gga_I_c[12] = {-0.0012170932750450447, 0.00559962844090478, -0.0253682622769702, 0.36712234112928854, -0.11636622422884102, -0.04014704257497215, -0.030403737273445586, 
  -0.01164221165527124, -0.02569053796469234, 0.05402566003613174, -0.11556358612142503, 0.0606093795504421};

  // J - scaled (for 2 elements per electron) GGA no erf contribution removed. Even temper generated.
  int al_gga_J_key = 104;
  int al_gga_J_ion_charge = 3;
  int al_gga_J_Ncoeff = 12;
  int al_gga_J_Nepe = 2;
  double al_gga_J_a[12] = {28.480358684358013, 10.0, 3.5111917342151315, 1.2328467394420664, 0.43287612810830595, 0.15199110829529341, 0.05336699231206311, 
  0.018738174228603847, 0.006579332246575677, 0.002310129700083161, 0.0008111308307896868, 0.0002848035868435803};
  double al_gga_J_c[12] = {-0.009128199555561878, 0.04199721329496242, -0.19026196702907328, 2.753417558444198, -0.8727466816781089, -0.30110281937231775, 
  -0.22802802952355705, -0.08731658762553707, -0.192679034313187, 0.4051924468949437, -0.866726890206337, 0.45457034185528755};

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
