#ifdef FIX_CLASS
FixStyle(fermi_nvt,FixFermi)
#else
#ifndef LMP_FIX_FERMI_H
#define LMP_FIX_FERMI_H

#include "fix.h"

namespace LAMMPS_NS {
class FixFermi : public Fix {
 public:
  FixFermi(class LAMMPS *, int, char **);
  virtual ~FixFermi();
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
	double dt;
	double dthalf;
	double dt4;
	double alpha_mass;
	double species_mass;
	double exp_cutoff;
	double t_target;
	double full_box_len;
	double alpha_half;
	double alpha_minus_quarter;
	double alpha_plus_quarter;
	double alpha_dot;
	double boltz;
	double mu;
	double initial_alpha;
	int old_nlocal;
	double** v_minus_quarter;
	double** v_plus_quarter;
	double** f_corrected;
    char *id_temp;

    class Compute *temperature,*pressure;

    int tcomputeflag;

};

}

#endif
#endif
