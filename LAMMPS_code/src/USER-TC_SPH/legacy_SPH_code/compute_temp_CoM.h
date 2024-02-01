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

#ifdef COMPUTE_CLASS

ComputeStyle(tem_CoM,ComputeTempCoM)

#else

#ifndef LMP_COMPUTE_TEMP_COM_H
#define LMP_COMPUTE_TEMP_COM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeTempCoM : public Compute {
 public:
  ComputeTempCoM(class LAMMPS *, int, char **);
  virtual ~ComputeTempCoM();
  void init() {}
  void setup();
  virtual double compute_scalar();
  virtual double compute_scalar_CoM(double N_epe, int N_ele, int tag_ele_start);
  virtual void compute_vector();

 protected:
  double tfactor;

  virtual void dof_compute();

//   double *vx_CoM;                  // CoM velocities (note not using atom parameter 'vx_COM')
//   double *vy_CoM;
//   double *vz_CoM;
  int N_ele_compute_temp;
  int size_vector_CoM_vel;
  // int CoM_vel_flag;
  int nmax;                        // required for communication of CoM velocity
  
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Temperature compute degrees of freedom < 0

This should not happen if you are calculating the temperature
on a valid set of atoms.

*/
