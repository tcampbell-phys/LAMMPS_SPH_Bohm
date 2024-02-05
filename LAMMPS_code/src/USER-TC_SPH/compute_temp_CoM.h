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

ComputeStyle(temp_CoM,ComputeTempCentreMass)

#else

#ifndef LMP_COMPUTE_TEMP_CENTRE_MASS_H
#define LMP_COMPUTE_TEMP_CENTRE_MASS_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeTempCentreMass : public Compute {
 public:
  ComputeTempCentreMass(class LAMMPS *, int, char **);
  virtual ~ComputeTempCentreMass();
  void init() {}
  void setup();
  virtual double compute_scalar();
  virtual double compute_scalar_CoM(double N_epe, double particle_mass, int N_ele, int tag_ele_start);
  virtual void compute_vector();

 protected:
  double tfactor;

  virtual void dof_compute();
  double *CoM_vel,*CoM_vel_all;
  int allocated;
  int unallocated;

  virtual void allocate();
  void deallocate();

  int N_ele_compute_temp;

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
