/* -*- c++ -*- ----------------------------------------------------------
   Gaussian decomposition parameters.

   For use with pair_coul_long_SPH within Bohm SPH framework.
   Thomas Campbell (Oxford)
------------------------------------------------------------------------- */
#ifndef LMP_PSEUDO_H
#define LMP_PSEUDO_H

namespace LAMMPS_NS {

class Pseudo {
 protected:
  // al.lda.lps BLPS

  // system A:
  // - N_ele = 3000
  // - L = 38.748135232811201
  // - sigma = 1.05
  char al_lda_A_name[] ='al_lda_A';
  int al_lda_A_N_ele = 3000;
  double al_lda_A_L = 38.748135232811201;
  double al_lda_A_sigma = 1.05;
  double al_lda_A_tfwhm = 13.285790472304193;
  int al_lda_A_Ncoeff = 15;
  static double *al_lda_A_c;
  static double *al_lda_A_a;
};

}

#endif