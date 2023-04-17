/* -*- c++ -*- ----------------------------------------------------------
Bespoke pair_style to compute the BohmGlobal force felt by electrons, according
to a Gaussian wave packet scheme.

Thomas Campbell
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(bohm_global_ions_full_ham,PairBohmGlobalIonsFullHam)

#else

#ifndef LMP_PAIR_BOHM_GLOBAL_IONS_FULL_HAM_H
#define LMP_PAIR_BOHM_GLOBAL_IONS_FULL_HAM_H

#include "pair.h"

namespace LAMMPS_NS {

class PairBohmGlobalIonsFullHam : public Pair {
 public:
  PairBohmGlobalIonsFullHam(class LAMMPS *);
  virtual ~PairBohmGlobalIonsFullHam();

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
  void *extract(const char *, int &);
 

protected:
  int nmax;                         // allocated size of per-atom arrays

  double e_sigma;
  double e_sigma_2;
  double e_sigma_3;
  double e_sigma_4;
  double ion_sigma;
  double ion_sigma_2;
  double ion_sigma_3;
  double ion_sigma_4;

  double f_prefactor;

  double e_packet_zeta_sigma2;
  double ion_packet_zeta_sigma2;

  double e_packet_pre;
  double ion_packet_pre;

  //hard code in electron mass for now

  double e_mass = 0.0005485799;

  double cut_global;
  double sigma_global;
  double ion_sigma_global;
  double ion_core_electrons;
  double gamma_factor;
  double hplanck;
  double hbar;
  double **cut;

  int *numforce;

  // per-atom arrays

  double *sum_packets;
    
  double *first_order_x_term;
  double *first_order_y_term;
  double *first_order_z_term;

  double *theta_x;
  double *theta_y;
  double *theta_z;

  double *first_order_x_full_ham;
  double *first_order_y_full_ham;
  double *first_order_z_full_ham;

  double *second_order_x_full_ham;
  double *second_order_y_full_ham;
  double *second_order_z_full_ham;

  double *third_order_x_full_ham;
  double *third_order_y_full_ham;
  double *third_order_z_full_ham;

  double *psi;

  double *rho_x;
  double *rho_y;
  double *rho_z;
  
  double *rho_x1;
  double *rho_x2;
  double *rho_x3;
  
  double *rho_y1;
  double *rho_y3;
  
  double *rho_z1;

  int commflag;


  virtual void allocate();
};

}

#endif
#endif