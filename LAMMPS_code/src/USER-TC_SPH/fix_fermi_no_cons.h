#ifdef FIX_CLASS
FixStyle(fermi_nvt_no_cons,FixFermiNoCons)
#else
#ifndef LMP_FIX_FERMI_NOCONS_H
#define LMP_FIX_FERMI_NOCONS_H

#include "fix.h"

namespace LAMMPS_NS {
class FixFermiNoCons : public Fix {
 public:
  FixFermiNoCons(class LAMMPS *, int, char **);
  int setmask();
  void init();
  void setup(int);
  void initial_integrate(int);
  void final_integrate();
  void write_restart(FILE *);
  virtual int pack_restart_data(double *); // pack restart data
  virtual void restart(char *);
  //double compute_scalar();

 protected:
	double dtf;
	double dtv,dtv_2;
	double alpha_mass;
	double species_mass;
	double exp_cutoff;
	double t_target;
	double alpha_half;
	double alpha_minus_quarter;
	double alpha_plus_quarter;
	double alpha_dot;
	double mu;
	double** v_minus_quarter;
	double** v_plus_quarter;
  char *id_temp;

  class Compute *temperature,*pressure;
  int tcomputeflag;

};

}

#endif
#endif
