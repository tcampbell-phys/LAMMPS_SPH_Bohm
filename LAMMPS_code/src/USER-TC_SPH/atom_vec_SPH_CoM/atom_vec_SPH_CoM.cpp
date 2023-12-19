/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   Edited by Thomas Campbell for use with dynamic width_SPHs in Bohm SPH scheme.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "atom_vec_SPH_CoM.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "memory.h"
#include "error.h"
#include "utils.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

AtomVecSPHCoM::AtomVecSPHCoM(LAMMPS *lmp) : AtomVec(lmp)
{
  molecular = 0;
  mass_type = 1;

  comm_x_only = 0;
  comm_f_only = 1;
  size_forward = 13;
  size_reverse = 3;
  size_border = 17;
  size_velocity = 3;
  size_data_atom = 16;
  size_data_vel = 4;
  xcol_data = 14;

  atom->q_flag = 1;
  atom->TC_SPH_flag = 1;
  atom->TC_X_CoM_SPH_flag = 1;
}

/* ----------------------------------------------------------------------
   grow atom arrays
   n = 0 grows arrays by a chunk
   n > 0 allocates arrays to size n
------------------------------------------------------------------------- */

void AtomVecSPHCoM::grow(int n)
{
  if (n == 0) grow_nmax();
  else nmax = n;
  atom->nmax = nmax;
  if (nmax < 0 || nmax > MAXSMALLINT)
    error->one(FLERR,"Per-processor system is too big");

  tag = memory->grow(atom->tag,nmax,"atom:tag");
  type = memory->grow(atom->type,nmax,"atom:type");
  mask = memory->grow(atom->mask,nmax,"atom:mask");
  image = memory->grow(atom->image,nmax,"atom:image");
  x = memory->grow(atom->x,nmax,3,"atom:x");
  v = memory->grow(atom->v,nmax,3,"atom:v");
  f = memory->grow(atom->f,nmax*comm->nthreads,3,"atom:f");

  //SPH density, width_SPH and omega_SPH terms:
  rho_SPH = memory->grow(atom->rho_SPH,nmax,"atom:rho_SPH");
  dx_rho_SPH = memory->grow(atom->dx_rho_SPH,nmax,"atom:dx_rho_SPH");
  dy_rho_SPH = memory->grow(atom->dy_rho_SPH,nmax,"atom:dy_rho_SPH");
  dz_rho_SPH = memory->grow(atom->dz_rho_SPH,nmax,"atom:dz_rho_SPH");
  width_SPH = memory->grow(atom->width_SPH,nmax,"atom:width_SPH");
  omega_SPH = memory->grow(atom->omega_SPH,nmax,"atom:omega_SPH");
  u_SPH = memory->grow(atom->u_SPH,nmax,"atom:u_SPH");

  //Bohm SPH CoM coords
  x_COM = memory->grow(atom->x_COM,nmax,"atom:x_COM");
  y_COM = memory->grow(atom->y_COM,nmax,"atom:y_COM");
  z_COM = memory->grow(atom->z_COM,nmax,"atom:z_COM");

  q = memory->grow(atom->q,nmax,"atom:q");

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->grow_arrays(nmax);
}

/* ----------------------------------------------------------------------
   reset local array ptrs
------------------------------------------------------------------------- */

void AtomVecSPHCoM::grow_reset()
{
  tag = atom->tag; type = atom->type;
  mask = atom->mask; image = atom->image;
  x = atom->x; v = atom->v; f = atom->f;
  q = atom->q;
  //SPH density, width_SPH and omega_SPH terms:
  rho_SPH = atom->rho_SPH; width_SPH = atom->width_SPH; omega_SPH = atom->omega_SPH; u_SPH = atom->u_SPH;
  dx_rho_SPH = atom->dx_rho_SPH; dy_rho_SPH = atom->dy_rho_SPH; dz_rho_SPH = atom->dz_rho_SPH;
  //Bohm SPH CoM coords:
  x_COM = atom->x_COM; y_COM = atom->y_COM; z_COM = atom->z_COM; 
}

/* ----------------------------------------------------------------------
   copy atom I info to atom J
------------------------------------------------------------------------- */

void AtomVecSPHCoM::copy(int i, int j, int delflag)
{
  tag[j] = tag[i];
  type[j] = type[i];
  mask[j] = mask[i];
  image[j] = image[i];
  x[j][0] = x[i][0];
  x[j][1] = x[i][1];
  x[j][2] = x[i][2];
  v[j][0] = v[i][0];
  v[j][1] = v[i][1];
  v[j][2] = v[i][2];

  q[j] = q[i];

  //SPH density, width_SPH and omega_SPH terms:
  rho_SPH[j] = rho_SPH[i];
  dx_rho_SPH[j] = dx_rho_SPH[i];
  dy_rho_SPH[j] = dy_rho_SPH[i];
  dz_rho_SPH[j] = dz_rho_SPH[i];
  width_SPH[j] = width_SPH[i];
  omega_SPH[j] = omega_SPH[i];
  u_SPH[j] = u_SPH[i];

  //Bohm SPH CoM coords:
  x_COM[j] = x_COM[i];
  y_COM[j] = y_COM[i];
  z_COM[j] = z_COM[i];

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->copy_arrays(i,j,delflag);
}

/* ---------------------------------------------------------------------- */

int AtomVecSPHCoM::pack_comm(int n, int *list, double *buf,
                             int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = rho_SPH[j];
      buf[m++] = dx_rho_SPH[j];
      buf[m++] = dy_rho_SPH[j];
      buf[m++] = dz_rho_SPH[j];
      buf[m++] = width_SPH[j];
      buf[m++] = omega_SPH[j];
      // set u_SPH to 0 on ghost processors to ensure consistency with
      // newton on
      buf[m++] = 0.0;
      //Bohm SPH CoM coords:
      buf[m++] = x_COM[j];
      buf[m++] = y_COM[j];
      buf[m++] = z_COM[j];

    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0]*domain->xprd + pbc[5]*domain->xy + pbc[4]*domain->xz;
      dy = pbc[1]*domain->yprd + pbc[3]*domain->yz;
      dz = pbc[2]*domain->zprd;
    }
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0] + dx;
      buf[m++] = x[j][1] + dy;
      buf[m++] = x[j][2] + dz;
      buf[m++] = rho_SPH[j];
      buf[m++] = dx_rho_SPH[j];
      buf[m++] = dy_rho_SPH[j];
      buf[m++] = dz_rho_SPH[j];
      buf[m++] = width_SPH[j];
      buf[m++] = omega_SPH[j];
      // set u_SPH to 0 on ghost processors to ensure consistency with
      // newton on
      buf[m++] = 0.0;
      //Bohm SPH CoM coords:
      buf[m++] = x_COM[j] + dx;
      buf[m++] = y_COM[j] + dy;
      buf[m++] = z_COM[j] + dz;

    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecSPHCoM::pack_comm_vel(int n, int *list, double *buf,
                                 int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz,dvx,dvy,dvz;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = rho_SPH[j];
      buf[m++] = dx_rho_SPH[j];
      buf[m++] = dy_rho_SPH[j];
      buf[m++] = dz_rho_SPH[j];
      buf[m++] = width_SPH[j];
      buf[m++] = omega_SPH[j];
      // set u_SPH to 0 on ghost processors to ensure consistency with
      // newton on
      buf[m++] = 0.0;
      //Bohm SPH CoM coords:
      buf[m++] = x_COM[j];
      buf[m++] = y_COM[j];
      buf[m++] = z_COM[j];
      buf[m++] = v[j][0];
      buf[m++] = v[j][1];
      buf[m++] = v[j][2];
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0]*domain->xprd + pbc[5]*domain->xy + pbc[4]*domain->xz;
      dy = pbc[1]*domain->yprd + pbc[3]*domain->yz;
      dz = pbc[2]*domain->zprd;
    }
    if (!deform_vremap) {
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;
        buf[m++] = rho_SPH[j];
        buf[m++] = dx_rho_SPH[j];
        buf[m++] = dy_rho_SPH[j];
        buf[m++] = dz_rho_SPH[j];
        buf[m++] = width_SPH[j];
        buf[m++] = omega_SPH[j];
        // set u_SPH to 0 on ghost processors to ensure consistency with
        // newton on
        buf[m++] = 0.0;
        //Bohm SPH CoM coords:
        buf[m++] = x_COM[j] + dx;
        buf[m++] = y_COM[j] + dy;
        buf[m++] = z_COM[j] + dz;
        buf[m++] = v[j][0];
        buf[m++] = v[j][1];
        buf[m++] = v[j][2];

      }
    } else {
      dvx = pbc[0]*h_rate[0] + pbc[5]*h_rate[5] + pbc[4]*h_rate[4];
      dvy = pbc[1]*h_rate[1] + pbc[3]*h_rate[3];
      dvz = pbc[2]*h_rate[2];
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;
        buf[m++] = rho_SPH[j];
        buf[m++] = dx_rho_SPH[j];
        buf[m++] = dy_rho_SPH[j];
        buf[m++] = dz_rho_SPH[j];
        buf[m++] = width_SPH[j];
        buf[m++] = omega_SPH[j];
        // set u_SPH to 0 on ghost processors to ensure consistency with
        // newton on
        buf[m++] = 0.0;
        //Bohm SPH CoM coords:
        buf[m++] = x_COM[j] + dx;
        buf[m++] = y_COM[j] + dy;
        buf[m++] = z_COM[j] + dz;
        if (mask[i] & deform_groupbit) {
          buf[m++] = v[j][0] + dvx;
          buf[m++] = v[j][1] + dvy;
          buf[m++] = v[j][2] + dvz;
        } else {
          buf[m++] = v[j][0];
          buf[m++] = v[j][1];
          buf[m++] = v[j][2];
        }
      }
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecSPHCoM::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    rho_SPH[i] = buf[m++];
    dx_rho_SPH[i] = buf[m++];
    dy_rho_SPH[i] = buf[m++];
    dz_rho_SPH[i] = buf[m++];
    width_SPH[i] = buf[m++];
    omega_SPH[i] = buf[m++];
    u_SPH[i] = buf[m++];
    //Bohm SPH CoM coords:
    x_COM[i] = buf[m++];
    y_COM[i] = buf[m++];
    z_COM[i] = buf[m++];

  }
}

/* ---------------------------------------------------------------------- */

void AtomVecSPHCoM::unpack_comm_vel(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    rho_SPH[i] = buf[m++];
    dx_rho_SPH[i] = buf[m++];
    dy_rho_SPH[i] = buf[m++];
    dz_rho_SPH[i] = buf[m++];
    width_SPH[i] = buf[m++];
    omega_SPH[i] = buf[m++];
    u_SPH[i] = buf[m++];
    //Bohm SPH CoM coords:
    x_COM[i] = buf[m++];
    y_COM[i] = buf[m++];
    z_COM[i] = buf[m++];
    v[i][0] = buf[m++];
    v[i][1] = buf[m++];
    v[i][2] = buf[m++];

  }
}

/* ---------------------------------------------------------------------- */

int AtomVecSPHCoM::pack_reverse(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = f[i][0];
    buf[m++] = f[i][1];
    buf[m++] = f[i][2];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecSPHCoM::unpack_reverse(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    f[j][0] += buf[m++];
    f[j][1] += buf[m++];
    f[j][2] += buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecSPHCoM::pack_border(int n, int *list, double *buf,
                               int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = rho_SPH[j];
      buf[m++] = dx_rho_SPH[j];
      buf[m++] = dy_rho_SPH[j];
      buf[m++] = dz_rho_SPH[j];
      buf[m++] = width_SPH[j];
      buf[m++] = omega_SPH[j];
      // set u_SPH to 0 on ghost processors to ensure consistency with
      // newton on
      buf[m++] = 0.0;
      //Bohm SPH CoM coords:
      buf[m++] = x_COM[j];
      buf[m++] = y_COM[j];
      buf[m++] = z_COM[j];
      buf[m++] = ubuf(tag[j]).d;
      buf[m++] = ubuf(type[j]).d;
      buf[m++] = ubuf(mask[j]).d;
      buf[m++] = q[j];
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0];
      dy = pbc[1];
      dz = pbc[2];
    }
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0] + dx;
      buf[m++] = x[j][1] + dy;
      buf[m++] = x[j][2] + dz;
      buf[m++] = rho_SPH[j];
      buf[m++] = dx_rho_SPH[j];
      buf[m++] = dy_rho_SPH[j];
      buf[m++] = dz_rho_SPH[j];
      buf[m++] = width_SPH[j];
      buf[m++] = omega_SPH[j];
      // set u_SPH to 0 on ghost processors to ensure consistency with
      // newton on
      buf[m++] = 0.0;
      //Bohm SPH CoM coords:
      buf[m++] = x_COM[j]+ dx;
      buf[m++] = y_COM[j]+ dy;
      buf[m++] = z_COM[j]+ dz;

      buf[m++] = ubuf(tag[j]).d;
      buf[m++] = ubuf(type[j]).d;
      buf[m++] = ubuf(mask[j]).d;
      buf[m++] = q[j];
    }
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->pack_border(n,list,&buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecSPHCoM::pack_border_vel(int n, int *list, double *buf,
                                   int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz,dvx,dvy,dvz;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = rho_SPH[j];
      buf[m++] = dx_rho_SPH[j];
      buf[m++] = dy_rho_SPH[j];
      buf[m++] = dz_rho_SPH[j];
      buf[m++] = width_SPH[j];
      buf[m++] = omega_SPH[j];
      // set u_SPH to 0 on ghost processors to ensure consistency with
      // newton on
      buf[m++] = 0.0;
      //Bohm SPH CoM coords:
      buf[m++] = x_COM[j];
      buf[m++] = y_COM[j];
      buf[m++] = z_COM[j];
      buf[m++] = ubuf(tag[j]).d;
      buf[m++] = ubuf(type[j]).d;
      buf[m++] = ubuf(mask[j]).d;
      buf[m++] = q[j];
      buf[m++] = v[j][0];
      buf[m++] = v[j][1];
      buf[m++] = v[j][2];
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0];
      dy = pbc[1];
      dz = pbc[2];
    }
    if (!deform_vremap) {
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;
        buf[m++] = rho_SPH[j];
        buf[m++] = dx_rho_SPH[j];
        buf[m++] = dy_rho_SPH[j];
        buf[m++] = dz_rho_SPH[j];
        buf[m++] = width_SPH[j];
        buf[m++] = omega_SPH[j];
        // set u_SPH to 0 on ghost processors to ensure consistency with
        // newton on
        buf[m++] = 0.0;
        //Bohm SPH CoM coords:
        buf[m++] = x_COM[j]+ dx;
        buf[m++] = y_COM[j]+ dy;
        buf[m++] = z_COM[j]+ dz;
        buf[m++] = ubuf(tag[j]).d;
        buf[m++] = ubuf(type[j]).d;
        buf[m++] = ubuf(mask[j]).d;
        buf[m++] = q[j];
        buf[m++] = v[j][0];
        buf[m++] = v[j][1];
        buf[m++] = v[j][2];
      }
    } else {
      dvx = pbc[0]*h_rate[0] + pbc[5]*h_rate[5] + pbc[4]*h_rate[4];
      dvy = pbc[1]*h_rate[1] + pbc[3]*h_rate[3];
      dvz = pbc[2]*h_rate[2];
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;
        buf[m++] = rho_SPH[j];
        buf[m++] = dx_rho_SPH[j];
        buf[m++] = dy_rho_SPH[j];
        buf[m++] = dz_rho_SPH[j];
        buf[m++] = width_SPH[j];
        buf[m++] = omega_SPH[j];
        // set u_SPH to 0 on ghost processors to ensure consistency with
        // newton on
        buf[m++] = 0.0;
        //Bohm SPH CoM coords:
        buf[m++] = x_COM[j]+ dx;
        buf[m++] = y_COM[j]+ dy;
        buf[m++] = z_COM[j]+ dz;
        buf[m++] = ubuf(tag[j]).d;
        buf[m++] = ubuf(type[j]).d;
        buf[m++] = ubuf(mask[j]).d;
        buf[m++] = q[j];
        if (mask[i] & deform_groupbit) {
          buf[m++] = v[j][0] + dvx;
          buf[m++] = v[j][1] + dvy;
          buf[m++] = v[j][2] + dvz;
        } else {
          buf[m++] = v[j][0];
          buf[m++] = v[j][1];
          buf[m++] = v[j][2];
        }
      }
    }
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->pack_border(n,list,&buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecSPHCoM::pack_border_hybrid(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = q[j];
    buf[m++] = rho_SPH[j];
    buf[m++] = dx_rho_SPH[j];
    buf[m++] = dy_rho_SPH[j];
    buf[m++] = dz_rho_SPH[j];
    buf[m++] = width_SPH[j];
    buf[m++] = omega_SPH[j];
    // set u_SPH to 0 on ghost processors to ensure consistency with
    // newton on
    buf[m++] = 0.0;
    //Bohm SPH CoM coords:
    buf[m++] = x_COM[j];
    buf[m++] = y_COM[j];
    buf[m++] = z_COM[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecSPHCoM::unpack_border(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (i == nmax) grow(0);
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    rho_SPH[i] = buf[m++];
    dx_rho_SPH[i] = buf[m++];
    dy_rho_SPH[i] = buf[m++];
    dz_rho_SPH[i] = buf[m++];
    width_SPH[i] = buf[m++];
    omega_SPH[i] = buf[m++];
    u_SPH[i] = buf[m++];
    //Bohm SPH CoM coords:
    x_COM[i] = buf[m++];
    y_COM[i] = buf[m++];
    z_COM[i] = buf[m++];
    tag[i] = (tagint) ubuf(buf[m++]).i;
    type[i] = (int) ubuf(buf[m++]).i;
    mask[i] = (int) ubuf(buf[m++]).i;
    q[i] = buf[m++];
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->
        unpack_border(n,first,&buf[m]);
}

/* ---------------------------------------------------------------------- */

void AtomVecSPHCoM::unpack_border_vel(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (i == nmax) grow(0);
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    rho_SPH[i] = buf[m++];
    dx_rho_SPH[i] = buf[m++];
    dy_rho_SPH[i] = buf[m++];
    dz_rho_SPH[i] = buf[m++];
    width_SPH[i] = buf[m++];
    omega_SPH[i] = buf[m++];
    u_SPH[i] = buf[m++];
    //Bohm SPH CoM coords:
    x_COM[i] = buf[m++];
    y_COM[i] = buf[m++];
    z_COM[i] = buf[m++];
    tag[i] = (tagint) ubuf(buf[m++]).i;
    type[i] = (int) ubuf(buf[m++]).i;
    mask[i] = (int) ubuf(buf[m++]).i;
    q[i] = buf[m++];
    v[i][0] = buf[m++];
    v[i][1] = buf[m++];
    v[i][2] = buf[m++];
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->
        unpack_border(n,first,&buf[m]);
}

/* ---------------------------------------------------------------------- */

int AtomVecSPHCoM::unpack_border_hybrid(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++)
    q[i] = buf[m++];
    rho_SPH[i] = buf[m++];
    dx_rho_SPH[i] = buf[m++];
    dy_rho_SPH[i] = buf[m++];
    dz_rho_SPH[i] = buf[m++];
    width_SPH[i] = buf[m++];
    omega_SPH[i] = buf[m++];
    u_SPH[i] = buf[m++];
    //Bohm SPH CoM coords:
    x_COM[i] = buf[m++];
    y_COM[i] = buf[m++];
    z_COM[i] = buf[m++];
  return m;
}

/* ----------------------------------------------------------------------
   pack data for atom I for sending to another proc
   xyz must be 1st 3 values, so comm::exchange() can test on them
------------------------------------------------------------------------- */

int AtomVecSPHCoM::pack_exchange(int i, double *buf)
{
  int m = 1;
  buf[m++] = x[i][0];
  buf[m++] = x[i][1];
  buf[m++] = x[i][2];
  buf[m++] = v[i][0];
  buf[m++] = v[i][1];
  buf[m++] = v[i][2];
  buf[m++] = ubuf(tag[i]).d;
  buf[m++] = ubuf(type[i]).d;
  buf[m++] = ubuf(mask[i]).d;
  buf[m++] = ubuf(image[i]).d;

  buf[m++] = q[i];
  buf[m++] = rho_SPH[i];
  buf[m++] = dx_rho_SPH[i];
  buf[m++] = dy_rho_SPH[i];
  buf[m++] = dz_rho_SPH[i];
  buf[m++] = width_SPH[i];
  buf[m++] = omega_SPH[i];
  buf[m++] = u_SPH[i];
  //Bohm SPH CoM coords:
  buf[m++] = x_COM[i];
  buf[m++] = y_COM[i];
  buf[m++] = z_COM[i];
  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      m += modify->fix[atom->extra_grow[iextra]]->pack_exchange(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecSPHCoM::unpack_exchange(double *buf)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  int m = 1;
  x[nlocal][0] = buf[m++];
  x[nlocal][1] = buf[m++];
  x[nlocal][2] = buf[m++];
  v[nlocal][0] = buf[m++];
  v[nlocal][1] = buf[m++];
  v[nlocal][2] = buf[m++];
  tag[nlocal] = (tagint) ubuf(buf[m++]).i;
  type[nlocal] = (int) ubuf(buf[m++]).i;
  mask[nlocal] = (int) ubuf(buf[m++]).i;
  image[nlocal] = (imageint) ubuf(buf[m++]).i;

  q[nlocal] = buf[m++];
  rho_SPH[nlocal] = buf[m++];
  dx_rho_SPH[nlocal] = buf[m++];
  dy_rho_SPH[nlocal] = buf[m++];
  dz_rho_SPH[nlocal] = buf[m++];
  width_SPH[nlocal] = buf[m++];
  omega_SPH[nlocal] = buf[m++];
  u_SPH[nlocal] = buf[m++];

  //Bohm SPH CoM coords:
  x_COM[nlocal] = buf[m++];
  y_COM[nlocal] = buf[m++];
  z_COM[nlocal] = buf[m++];


  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      m += modify->fix[atom->extra_grow[iextra]]->
        unpack_exchange(nlocal,&buf[m]);

  atom->nlocal++;
  return m;
}

/* ----------------------------------------------------------------------
   size of restart data for all atoms owned by this proc
   include extra data stored by fixes
------------------------------------------------------------------------- */

int AtomVecSPHCoM::size_restart()
{
  int i;

  int nlocal = atom->nlocal;

  int n = 22 * nlocal;

  if (atom->nextra_restart)
    for (int iextra = 0; iextra < atom->nextra_restart; iextra++)
      for (i = 0; i < nlocal; i++)
        n += modify->fix[atom->extra_restart[iextra]]->size_restart(i);

  return n;
}

/* ----------------------------------------------------------------------
   pack atom I's data for restart file including extra quantities
   xyz must be 1st 3 values, so that read_restart can test on them
   molecular types may be negative, but write as positive
------------------------------------------------------------------------- */

int AtomVecSPHCoM::pack_restart(int i, double *buf)
{
  int m = 1;
  buf[m++] = x[i][0];
  buf[m++] = x[i][1];
  buf[m++] = x[i][2];
  buf[m++] = rho_SPH[i];
  buf[m++] = dx_rho_SPH[i];
  buf[m++] = dy_rho_SPH[i];
  buf[m++] = dz_rho_SPH[i];
  buf[m++] = width_SPH[i];
  buf[m++] = omega_SPH[i];
  buf[m++] = u_SPH[i];
  //Bohm SPH CoM coords:
  buf[m++] = x_COM[i];
  buf[m++] = y_COM[i];
  buf[m++] = z_COM[i];
  buf[m++] = ubuf(tag[i]).d;
  buf[m++] = ubuf(type[i]).d;
  buf[m++] = ubuf(mask[i]).d;
  buf[m++] = ubuf(image[i]).d;
  buf[m++] = v[i][0];
  buf[m++] = v[i][1];
  buf[m++] = v[i][2];

  buf[m++] = q[i];

  if (atom->nextra_restart)
    for (int iextra = 0; iextra < atom->nextra_restart; iextra++)
      m += modify->fix[atom->extra_restart[iextra]]->pack_restart(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ----------------------------------------------------------------------
   unpack data for one atom from restart file including extra quantities
------------------------------------------------------------------------- */

int AtomVecSPHCoM::unpack_restart(double *buf)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) {
    grow(0);
    if (atom->nextra_store)
      memory->grow(atom->extra,nmax,atom->nextra_store,"atom:extra");
  }

  int m = 1;
  x[nlocal][0] = buf[m++];
  x[nlocal][1] = buf[m++];
  x[nlocal][2] = buf[m++];
  rho_SPH[nlocal] = buf[m++];
  dx_rho_SPH[nlocal] = buf[m++];
  dy_rho_SPH[nlocal] = buf[m++];
  dz_rho_SPH[nlocal] = buf[m++];
  width_SPH[nlocal] = buf[m++];
  omega_SPH[nlocal] = buf[m++];
  u_SPH[nlocal] = buf[m++];
  //Bohm SPH CoM coords:
  x_COM[nlocal] = buf[m++];
  y_COM[nlocal] = buf[m++];
  z_COM[nlocal] = buf[m++];
  tag[nlocal] = (tagint) ubuf(buf[m++]).i;
  type[nlocal] = (int) ubuf(buf[m++]).i;
  mask[nlocal] = (int) ubuf(buf[m++]).i;
  image[nlocal] = (imageint) ubuf(buf[m++]).i;
  v[nlocal][0] = buf[m++];
  v[nlocal][1] = buf[m++];
  v[nlocal][2] = buf[m++];

  q[nlocal] = buf[m++];
  

  double **extra = atom->extra;
  if (atom->nextra_store) {
    int size = static_cast<int> (buf[0]) - m;
    for (int i = 0; i < size; i++) extra[nlocal][i] = buf[m++];
  }

  atom->nlocal++;
  return m;
}

/* ----------------------------------------------------------------------
   create one atom of itype at coord
   set other values to defaults
------------------------------------------------------------------------- */

void AtomVecSPHCoM::create_atom(int itype, double *coord)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  tag[nlocal] = 0;
  type[nlocal] = itype;
  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];

  rho_SPH[nlocal] = 0.0;
  dx_rho_SPH[nlocal] = 0.0;
  dy_rho_SPH[nlocal] = 0.0;
  dz_rho_SPH[nlocal] = 0.0;
  width_SPH[nlocal] = 0.0;
  omega_SPH[nlocal] = 0.0;
  u_SPH[nlocal] = 0.0;
  //Bohm SPH CoM coords:
  x_COM[nlocal] = 0.0;
  y_COM[nlocal] = 0.0;
  z_COM[nlocal] = 0.0;

  mask[nlocal] = 1;
  image[nlocal] = ((imageint) IMGMAX << IMG2BITS) |
    ((imageint) IMGMAX << IMGBITS) | IMGMAX;
  v[nlocal][0] = 0.0;
  v[nlocal][1] = 0.0;
  v[nlocal][2] = 0.0;

  q[nlocal] = 0.0;

  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   unpack one line from Atoms section of data file
   initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecSPHCoM::data_atom(double *coord, imageint imagetmp, char **values)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  tag[nlocal] = utils::tnumeric(FLERR,values[0],true,lmp);
  type[nlocal] = utils::inumeric(FLERR,values[1],true,lmp);
  if (type[nlocal] <= 0 || type[nlocal] > atom->ntypes)
    error->one(FLERR,"Invalid atom type in Atoms section of data file");

  q[nlocal] = utils::numeric(FLERR,values[2],true,lmp);
  rho_SPH[nlocal] = utils::numeric(FLERR,values[3],true,lmp);
  dx_rho_SPH[nlocal] = utils::numeric(FLERR,values[4],true,lmp);
  dy_rho_SPH[nlocal] = utils::numeric(FLERR,values[5],true,lmp);
  dz_rho_SPH[nlocal] = utils::numeric(FLERR,values[6],true,lmp);
  width_SPH[nlocal] = utils::numeric(FLERR,values[7],true,lmp);
  omega_SPH[nlocal] = utils::numeric(FLERR,values[8],true,lmp);
  u_SPH[nlocal] = utils::numeric(FLERR,values[9],true,lmp);

  //Bohm SPH CoM coords:
  x_COM[nlocal] = utils::numeric(FLERR,values[10],true,lmp);
  y_COM[nlocal] = utils::numeric(FLERR,values[11],true,lmp);
  z_COM[nlocal] = utils::numeric(FLERR,values[12],true,lmp);

  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];

  image[nlocal] = imagetmp;

  mask[nlocal] = 1;
  v[nlocal][0] = 0.0;
  v[nlocal][1] = 0.0;
  v[nlocal][2] = 0.0;

  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   unpack hybrid quantities from one line in Atoms section of data file
   initialize other atom quantities for this sub-style
------------------------------------------------------------------------- */

int AtomVecSPHCoM::data_atom_hybrid(int nlocal, char **values)
{
  q[nlocal] = utils::numeric(FLERR,values[0],true,lmp);
  rho_SPH[nlocal] = 0.0;
  dx_rho_SPH[nlocal] = 0.0;
  dy_rho_SPH[nlocal] = 0.0;
  dz_rho_SPH[nlocal] = 0.0;
  width_SPH[nlocal] = 0.0;
  omega_SPH[nlocal] = 0.0;
  u_SPH[nlocal] = 0.0;
  //Bohm SPH CoM coords:
  x_COM[nlocal] = 0.0;
  y_COM[nlocal] = 0.0;
  z_COM[nlocal] = 0.0;

  return 11;
}

/* ----------------------------------------------------------------------
   pack atom info for data file including 3 image flags
------------------------------------------------------------------------- */

void AtomVecSPHCoM::pack_data(double **buf)
{
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    buf[i][0] = ubuf(tag[i]).d;
    buf[i][1] = ubuf(type[i]).d;
    buf[i][2] = q[i];
    buf[i][3] = rho_SPH[i];
    buf[i][4] = dx_rho_SPH[i];
    buf[i][5] = dy_rho_SPH[i];
    buf[i][6] = dz_rho_SPH[i];
    buf[i][7] = width_SPH[i];
    buf[i][8] = omega_SPH[i];
    buf[i][9] = u_SPH[i];
    //Bohm SPH CoM coords:
    buf[i][10] = x_COM[i];
    buf[i][11] = y_COM[i];
    buf[i][12] = z_COM[i];
    buf[i][13] = x[i][0];
    buf[i][14] = x[i][1];
    buf[i][15] = x[i][2];
    buf[i][16] = ubuf((image[i] & IMGMASK) - IMGMAX).d;
    buf[i][17] = ubuf((image[i] >> IMGBITS & IMGMASK) - IMGMAX).d;
    buf[i][18] = ubuf((image[i] >> IMG2BITS) - IMGMAX).d;
  }
}

/* ----------------------------------------------------------------------
   pack hybrid atom info for data file
------------------------------------------------------------------------- */

int AtomVecSPHCoM::pack_data_hybrid(int i, double *buf)
{
  buf[0] = q[i];
  buf[1] = rho_SPH[i];
  buf[2] = dx_rho_SPH[i];
  buf[3] = dy_rho_SPH[i];
  buf[4] = dz_rho_SPH[i];
  buf[5] = width_SPH[i];
  buf[6] = omega_SPH[i];
  buf[7] = u_SPH[i];
  //Bohm SPH CoM coords:
  buf[8] = x_COM[i];
  buf[9] = y_COM[i];
  buf[10] = z_COM[i];
  return 11;
}

/* ----------------------------------------------------------------------
   write atom info to data file including 3 image flags
------------------------------------------------------------------------- */

void AtomVecSPHCoM::write_data(FILE *fp, int n, double **buf)
{
  for (int i = 0; i < n; i++)
    fprintf(fp,TAGINT_FORMAT " %d %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %d %d %d\n",
            (tagint) ubuf(buf[i][0]).i,(int) ubuf(buf[i][1]).i,
            buf[i][2],buf[i][3],buf[i][4],buf[i][5],buf[i][6],buf[i][7],buf[i][8],buf[i][9],
            buf[i][10],buf[i][11],buf[i][12],buf[i][13],buf[i][14],buf[i][15],(int) ubuf(buf[i][16]).i,(int) ubuf(buf[i][17]).i,
            (int) ubuf(buf[i][18]).i);
}

/* ----------------------------------------------------------------------
   write hybrid atom info to data file
------------------------------------------------------------------------- */

int AtomVecSPHCoM::write_data_hybrid(FILE *fp, double *buf)
{
  fprintf(fp," %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e",buf[0],buf[1],buf[2],buf[3],buf[4],buf[5],buf[6],buf[7],buf[8],buf[9],buf[10]);
  return 11;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

bigint AtomVecSPHCoM::memory_usage()
{
  bigint bytes = 0;

  if (atom->memcheck("tag")) bytes += memory->usage(tag,nmax);
  if (atom->memcheck("type")) bytes += memory->usage(type,nmax);
  if (atom->memcheck("mask")) bytes += memory->usage(mask,nmax);
  if (atom->memcheck("image")) bytes += memory->usage(image,nmax);
  if (atom->memcheck("x")) bytes += memory->usage(x,nmax,3);
  if (atom->memcheck("v")) bytes += memory->usage(v,nmax,3);
  if (atom->memcheck("f")) bytes += memory->usage(f,nmax*comm->nthreads,3);

  if (atom->memcheck("q")) bytes += memory->usage(q,nmax);
  if (atom->memcheck("rho_SPH")) bytes += memory->usage(rho_SPH,nmax);
  if (atom->memcheck("dx_rho_SPH")) bytes += memory->usage(dx_rho_SPH,nmax);    
  if (atom->memcheck("dy_rho_SPH")) bytes += memory->usage(dy_rho_SPH,nmax);
  if (atom->memcheck("dz_rho_SPH")) bytes += memory->usage(dz_rho_SPH,nmax);
  if (atom->memcheck("width_SPH")) bytes += memory->usage(width_SPH,nmax);
  if (atom->memcheck("omega_SPH")) bytes += memory->usage(omega_SPH,nmax);
  if (atom->memcheck("u_SPH")) bytes += memory->usage(u_SPH,nmax);

  //Bohm SPH CoM coords:
  if (atom->memcheck("x_COM")) bytes += memory->usage(x_COM,nmax);
  if (atom->memcheck("y_COM")) bytes += memory->usage(y_COM,nmax);
  if (atom->memcheck("z_COM")) bytes += memory->usage(z_COM,nmax);

  return bytes;
}