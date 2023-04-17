#ifdef FIX_CLASS
FixStyle(fermi_nvt,FixFermi)
#else

class FixFermi : public Fix {
 public:
  FixFermi(class LAMMPS *, int, char **);
  virtual void init();
  virtual void initial_integrate();
  virtual void final_integrate();

 protected:
	double dt;
	double dthalf;
	double dt4;
	double alpha_mass;
	double alpha_half;
	double alpha_minus_quarter;
	double alpha_plus_quarter;
	double alpha_dot;
	double mu;
	double** v_minus_quarter;
	double** v_plus_quarter;