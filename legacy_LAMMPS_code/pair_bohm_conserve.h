/* -*- c++ -*- ----------------------------------------------------------
Bespoke pair_style to compute the BohmConserve force felt by electrons, according
to a Gaussian wave packet scheme.

Thomas Campbell
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(bohm_conserve,PairBohmConserve)

#else

#ifndef LMP_PAIR_BOHM_CONSERVE_H
#define LMP_PAIR_BOHM_CONSERVE_H

#include "pair.h"

namespace LAMMPS_NS {

class PairBohmConserve : public Pair {
 public:
  PairBohmConserve(class LAMMPS *);
  virtual ~PairBohmConserve();

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

  double cut_global;
  double hplanck;
  double hbar;
  double **cut;

  int *numforce;

  // per-atom arrays

  double *sum_packets;
    
  double *first_order_x_term;
  double *first_order_y_term;
  double *first_order_z_term;
  
  double *second_order_x_term;
  double *second_order_y_term;
  double *second_order_z_term;

  double *third_order_x_term;
  double *third_order_y_term;
  double *third_order_z_term;
  
  double *theta_x;
  double *theta_y;
  double *theta_z;

  double *phi;

  double *eta;

  double *rho_x;
  double *rho_y;
  double *rho_z;
  
  double *rho_x1;
  double *rho_x2;
  double *rho_x3;
  
  double *rho_y1;
  double *rho_y3;
  
  double *rho_z1;

  double *first_order_sder;
  double *second_order_sder;
  double *third_order_sder;

  int commflag;


  virtual void allocate();
};

}

#endif
#endif