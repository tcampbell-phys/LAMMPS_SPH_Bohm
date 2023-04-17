/* -*- c++ -*- ----------------------------------------------------------
Bespoke pair_style to compute the BohmGlobal force felt by electrons, according
to a Gaussian wave packet scheme.

Thomas Campbell
------------------------------------------------------------------------- */
#include "pair_bohm_global_ions.h"
#include <mpi.h>
#include <cmath>
#include <cstring>
#include <math.h>
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"
#include "update.h"
#include "utils.h"

using namespace LAMMPS_NS;

#define MAXLINE 1024

/* ---------------------------------------------------------------------- */

PairBohmGlobalIons::PairBohmGlobalIons(LAMMPS *lmp) : Pair(lmp) {
  nmax = 0;
  //fprintf(screen,"Variables in setup stage \n");
  manybody_flag = 1;

  sum_packets = NULL;
    
  first_order_x_term = NULL;
  first_order_y_term = NULL;
  first_order_z_term = NULL;
  
  theta_x = NULL;
  theta_y = NULL;
  theta_z = NULL;

  psi = NULL;

  rho_x = NULL;
  rho_y = NULL;
  rho_z = NULL;
  
  rho_x1 = NULL;
  rho_x2 = NULL;
  rho_x3 = NULL;
  
  rho_y1 = NULL;
  rho_y3 = NULL;
  
  rho_z1 = NULL;

  comm_forward = 17;
  comm_reverse = 17;
}

/* ---------------------------------------------------------------------- */

PairBohmGlobalIons::~PairBohmGlobalIons()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);

    memory->destroy(sum_packets);

    memory->destroy(first_order_x_term);
    memory->destroy(first_order_y_term);
    memory->destroy(first_order_z_term);
    
    memory->destroy(theta_x);
    memory->destroy(theta_y);
    memory->destroy(theta_z);

    memory->destroy(psi);

    memory->destroy(rho_x);
    memory->destroy(rho_y);
    memory->destroy(rho_z);
    
    memory->destroy(rho_x1);
    memory->destroy(rho_x2);
    memory->destroy(rho_x3);
    
    memory->destroy(rho_y1);
    memory->destroy(rho_y3);
    
    memory->destroy(rho_z1);
  }
}

/* ---------------------------------------------------------------------- */

void PairBohmGlobalIons::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz;
  double xtmp2,ytmp2,ztmp2,delx2,dely2,delz2;
  double delx_2,dely_2,delz_2;
  double delx_sigma2,dely_sigma2,delz_sigma2;
  double n_delx_sigma2,n_dely_sigma2,n_delz_sigma2;
  double sigma;
  double sigma_4,sigma_2;
  double phi_ij;
  double phi_ji;
  double omega_ij;
  double omega_ij2;
  double omega_ji;
  double sum_packets_i_2;
  double sum_packets_i_3;
  double rsq;
  double rsq2;
  double rsq_sigma2;
  double n_rsq_sigma2;
  double reduced_theta_x,reduced_theta_y,reduced_theta_z;
  double reduced_theta_x_i,reduced_theta_y_i,reduced_theta_z_i;
  double bohm_pot;
  double second_order_x_term_i,second_order_y_term_i,second_order_z_term_i;
  double third_order_x_term_i,third_order_y_term_i,third_order_z_term_i;
  double gamma_i;

  int *ilist,*jlist,*numneigh,**firstneigh;

  ev_init(eflag,vflag);

  if (atom->nmax > nmax) {
    memory->destroy(sum_packets);

    memory->destroy(first_order_x_term);
    memory->destroy(first_order_y_term);
    memory->destroy(first_order_z_term);
    
    memory->destroy(theta_x);
    memory->destroy(theta_y);
    memory->destroy(theta_z);

    memory->destroy(psi);

    memory->destroy(rho_x);
    memory->destroy(rho_y);
    memory->destroy(rho_z);
    
    memory->destroy(rho_x1);
    memory->destroy(rho_x2);
    memory->destroy(rho_x3);
    
    memory->destroy(rho_y1);
    memory->destroy(rho_y3);
    
    memory->destroy(rho_z1);

    nmax = atom->nmax;

    memory->create(sum_packets,nmax,"pair:sum_packets");

    memory->create(first_order_x_term,nmax,"pair:first_order_x_term");
    memory->create(first_order_y_term,nmax,"pair:first_order_y_term");
    memory->create(first_order_z_term,nmax,"pair:first_order_z_term");
    
    memory->create(theta_x,nmax,"pair:theta_x");
    memory->create(theta_y,nmax,"pair:theta_y");
    memory->create(theta_z,nmax,"pair:theta_z:");

    memory->create(psi,nmax,"pair:psi");

    memory->create(rho_x,nmax,"pair:rho_x");
    memory->create(rho_y,nmax,"pair:rho_y");
    memory->create(rho_z,nmax,"pair:rho_z");
    
    memory->create(rho_x1,nmax,"pair:rho_x1");
    memory->create(rho_x2,nmax,"pair:rho_x2");
    memory->create(rho_x3,nmax,"pair:rho_x3");
    
    memory->create(rho_y1,nmax,"pair:rho_y1");
    memory->create(rho_y3,nmax,"pair:rho_y3");
    
    memory->create(rho_z1,nmax,"pair:rho_z1");
  }

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;

  int *tag =atom->tag;

  double *mass = atom->mass;

  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int newton_pair = force->newton_pair;

  double cutsquared = 0.;

  double jloop_sigma;
  double jloop_sigma_2;
  double jloop_packet_pre;
  double jloop_packet_zeta_sigma2;

  double iloop_sigma;
  double iloop_sigma_2;
  double iloop_packet_pre;
  double iloop_packet_zeta_sigma2;

  double *sum_packets_d = atom->sum_packets_d;
  double *first_order_x_term_d = atom->first_order_x_term_d;
  double *first_order_y_term_d = atom->first_order_y_term_d;
  double *first_order_z_term_d = atom->first_order_z_term_d;
  double *phi_d = atom->phi_d;
  double *rho_x_d = atom->rho_x_d;
  double *rho_y_d = atom->rho_y_d;
  double *rho_z_d = atom->rho_z_d;
  double *theta_x_d = atom->theta_x_d;
  double *theta_y_d = atom->theta_y_d;
  double *theta_z_d = atom->theta_z_d;
  double *second_order_x_term_d = atom->second_order_x_term_d;
  double *second_order_y_term_d = atom->second_order_y_term_d;
  double *second_order_z_term_d = atom->second_order_z_term_d;
  double *third_order_x_term_d = atom->third_order_x_term_d;
  double *third_order_y_term_d = atom->third_order_y_term_d;
  double *third_order_z_term_d = atom->third_order_z_term_d;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  cutsquared = cut_global*cut_global;
  // zero out per-atom arrays

  if (newton_pair) {
    for (i = 0; i < nall; i++){
      sum_packets[i] = 0.0;
        
      first_order_x_term[i] = 0.0;
      first_order_y_term[i] = 0.0;
      first_order_z_term[i] = 0.0;
      
      theta_x[i] = 0.0;
      theta_y[i] = 0.0;
      theta_z[i] = 0.0;

      psi[i] = 0.0;

      rho_x[i] = 0.0;
      rho_y[i] = 0.0;
      rho_z[i] = 0.0;
      
      rho_x1[i] = 0.0;
      rho_x2[i] = 0.0;
      rho_x3[i] = 0.0;
      
      rho_y1[i] = 0.0;
      rho_y3[i] = 0.0;
      
      rho_z1[i] = 0.0;
    }
  } 
  else{
    for (i = 0; i < nlocal; i++){
      sum_packets[i] = 0.0;
        
      first_order_x_term[i] = 0.0;
      first_order_y_term[i] = 0.0;
      first_order_z_term[i] = 0.0;
      
      theta_x[i] = 0.0;
      theta_y[i] = 0.0;
      theta_z[i] = 0.0;

      psi[i] = 0.0;

      rho_x[i] = 0.0;
      rho_y[i] = 0.0;
      rho_z[i] = 0.0;
      
      rho_x1[i] = 0.0;
      rho_x2[i] = 0.0;
      rho_x3[i] = 0.0;
      
      rho_y1[i] = 0.0;
      rho_y3[i] = 0.0;
      
      rho_z1[i] = 0.0;
    }
  }

  // loop over my atoms

  /*NOTE - ions must be type 1 and electrons 2 or 3 in this implementation.*/

  /* Scheme for calculating the Bohm force with core electrons included explicitly as a Gaussian
  attached to each ion. */

  for (ii = 0; ii < inum; ii++) {

    i = ilist[ii];

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    itype = type[i];
    jlist = firstneigh[i];

    ////fprintf(screen,"itype = %d \n",itype);

    if (itype == 1){
      iloop_sigma = ion_sigma;
      ////fprintf(screen,"iloop_sigma = %f \n",iloop_sigma);
      iloop_sigma_2 = ion_sigma_2;
      iloop_packet_pre = ion_packet_pre;
      iloop_packet_zeta_sigma2 = ion_packet_zeta_sigma2;
    }

    if (itype == 2 || itype == 3){
      iloop_sigma = e_sigma;
      ////fprintf(screen,"iloop_sigma = %f \n",iloop_sigma);
      iloop_sigma_2 = e_sigma_2;
      iloop_packet_pre = e_packet_pre;
      iloop_packet_zeta_sigma2 = e_packet_zeta_sigma2;
    }

    jnum = numneigh[i];

     /* Self interaction terms */

    sum_packets[i] += iloop_packet_pre;


    /*N.B. 'psi' term below refers to the 'phi term from old Bohm implementations. Phi has been updated to it's meaning
    in 27/03/22 derivation.*/

    psi[i] += -3*iloop_packet_zeta_sigma2;
    rho_x1[i] += -iloop_packet_zeta_sigma2;
    rho_y1[i] += -iloop_packet_zeta_sigma2;
    rho_z1[i] += -iloop_packet_zeta_sigma2;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];

      rsq = delx*delx + dely*dely + delz*delz;

      jtype = type[j];

      
      ////fprintf(screen,"type[j] = %i \n",type[j]);
      if (rsq < cutsquared) {

        ////fprintf(screen,"jtype = %d \n",jtype);

        if (jtype == 1){
          jloop_sigma = ion_sigma;
          ////fprintf(screen,"jloop_sigma = %f \n",jloop_sigma);
          jloop_sigma_2 = ion_sigma_2;
          jloop_packet_pre = ion_packet_pre;
          jloop_packet_zeta_sigma2 = ion_packet_zeta_sigma2;
        }

        if (jtype == 2 || jtype == 3){
          jloop_sigma = e_sigma;
          ////fprintf(screen,"jloop_sigma = %f \n",jloop_sigma);
          jloop_sigma_2 = e_sigma_2;
          jloop_packet_pre = e_packet_pre;
          jloop_packet_zeta_sigma2 = e_packet_zeta_sigma2;
        }

        delx_2 = pow(delx,2);
        dely_2 = pow(dely,2);
        delz_2 = pow(delz,2);

        rsq_sigma2 = rsq/jloop_sigma_2;

        delx_sigma2 = delx/jloop_sigma_2;
        dely_sigma2 = dely/jloop_sigma_2;
        delz_sigma2 = delz/jloop_sigma_2;

        omega_ij = exp(-(rsq)/(2*jloop_sigma_2));

        phi_ij = jloop_packet_zeta_sigma2*omega_ij;

        sum_packets[i] += jloop_packet_pre*omega_ij;

        first_order_x_term[i] += phi_ij*(delx_sigma2)*(5 - rsq_sigma2);
        first_order_y_term[i] += phi_ij*(dely_sigma2)*(5 - rsq_sigma2);
        first_order_z_term[i] += phi_ij*(delz_sigma2)*(5 - rsq_sigma2);

        theta_x[i] += phi_ij*delx;
        theta_y[i] += phi_ij*dely;
        theta_z[i] += phi_ij*delz;

        psi[i] += phi_ij*(rsq_sigma2 - 3);

        rho_x1[i] += (-1 + delx_2/jloop_sigma_2)*(phi_ij);
        rho_x2[i] += (delx*dely)*(phi_ij)/jloop_sigma_2;
        rho_x3[i] += (delx*delz)*(phi_ij)/jloop_sigma_2;
        rho_y1[i] += (-1 + dely_2/jloop_sigma_2)*(phi_ij);
        rho_y3[i] += (dely*delz)*(phi_ij)/jloop_sigma_2;
        rho_z1[i] += (-1 + delz_2/jloop_sigma_2)*(phi_ij);

        if (newton_pair || j < nlocal) {

          omega_ji = exp(-(rsq)/(2*iloop_sigma_2));

          sum_packets[j] += iloop_packet_pre*omega_ji;

          phi_ji = iloop_packet_zeta_sigma2*omega_ji;

          /*'n_' prefix denotes newton quantities*/

          n_rsq_sigma2 = rsq/iloop_sigma_2;

          n_delx_sigma2 = delx/iloop_sigma_2;
          n_dely_sigma2 = dely/iloop_sigma_2;
          n_delz_sigma2 = delz/iloop_sigma_2;

          first_order_x_term[j] -= phi_ji*n_delx_sigma2*(5-n_rsq_sigma2);
          first_order_y_term[j] -= phi_ji*n_dely_sigma2*(5-n_rsq_sigma2);
          first_order_z_term[j] -= phi_ji*n_delz_sigma2*(5-n_rsq_sigma2);
          
          theta_x[j] -= phi_ji*delx;
          theta_y[j] -= phi_ji*dely;
          theta_z[j] -= phi_ji*delz;

          psi[j] += phi_ji*(n_rsq_sigma2-3);

          rho_x1[j] += (-1 + delx_2/iloop_sigma_2)*(phi_ji);
          rho_x2[j] += ((delx*dely)/iloop_sigma_2)*(phi_ji);
          rho_x3[j] += ((delx*delz)/iloop_sigma_2)*(phi_ji);
          rho_y1[j] += (-1 + dely_2/iloop_sigma_2)*(phi_ji);
          rho_y3[j] += ((dely*delz)/iloop_sigma_2)*(phi_ji);
          rho_z1[j] += (-1 + delz_2/iloop_sigma_2)*(phi_ji);

        }
      }
    }
  }
  
  commflag = 0;

  if (newton_pair) comm->reverse_comm_pair(this);

  comm->forward_comm_pair(this);

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    
    sum_packets_i_2 = pow(sum_packets[i],2);
    sum_packets_i_3 = pow(sum_packets[i],3);

    if (atom->forceterms_flag == 1){
      sum_packets_d[i] = sum_packets[i];
      first_order_x_term_d[i] = first_order_x_term[i];
      first_order_y_term_d[i] = first_order_y_term[i];
      first_order_z_term_d[i] = first_order_z_term[i];
    }

    first_order_x_term[i] *= (1/sum_packets[i]);
    first_order_y_term[i] *= (1/sum_packets[i]);
    first_order_z_term[i] *= (1/sum_packets[i]);

    rho_x[i] = theta_x[i]*rho_x1[i] + theta_y[i]*rho_x2[i] + theta_z[i]*rho_x3[i];
    
    rho_y[i] = theta_x[i]*rho_x2[i] + theta_y[i]*rho_y1[i] + theta_z[i]*rho_y3[i];

    rho_z[i] = theta_x[i]*rho_x3[i] + theta_y[i]*rho_y3[i] + theta_z[i]*rho_z1[i];

    second_order_x_term_i = psi[i]*theta_x[i] + rho_x[i];
    second_order_y_term_i = psi[i]*theta_y[i] + rho_y[i];
    second_order_z_term_i = psi[i]*theta_z[i] + rho_z[i];

    if (atom->forceterms_flag == 1){
      phi_d[i] = psi[i];

      rho_x_d[i] = rho_x[i];
      rho_y_d[i] = rho_y[i];
      rho_z_d[i] = rho_z[i];

      theta_x_d[i] = theta_x[i];
      theta_y_d[i] = theta_y[i];
      theta_z_d[i] = theta_z[i];

      second_order_x_term_d[i] = second_order_x_term_i;
      second_order_y_term_d[i] = second_order_y_term_i;
      second_order_z_term_d[i] = second_order_z_term_i;
    }
    
    second_order_x_term_i *= (1/(sum_packets_i_2));
    second_order_y_term_i *= (1/(sum_packets_i_2));
    second_order_z_term_i *= (1/(sum_packets_i_2));

    gamma_i = pow(theta_x[i],2) + pow(theta_y[i],2) + pow(theta_z[i],2);
    
    third_order_x_term_i = -theta_x[i]*gamma_i;
    third_order_y_term_i = -theta_y[i]*gamma_i;
    third_order_z_term_i = -theta_z[i]*gamma_i;

    if (atom->forceterms_flag == 1){
      third_order_x_term_d[i] = third_order_x_term_i;
      third_order_y_term_d[i] = third_order_y_term_i;
      third_order_z_term_d[i] = third_order_z_term_i;
    }
    
    third_order_x_term_i *= (1/(sum_packets_i_3));
    third_order_y_term_i *= (1/(sum_packets_i_3));
    third_order_z_term_i *= (1/(sum_packets_i_3));

    //fprintf(screen,"i = %i \n",i);

    //fprintf(screen,"first_order_x_term[i] = %f \n",first_order_x_term[i]);
    //fprintf(screen,"first_order_y_term[i] = %f \n",first_order_y_term[i]);
    //fprintf(screen,"first_order_z_term[i] = %f \n",first_order_z_term[i]);

    //fprintf(screen,"second_order_x_term_i = %f \n",second_order_x_term_i);
    //fprintf(screen,"second_order_y_term_i = %f \n",second_order_y_term_i);
    //fprintf(screen,"second_order_z_term_i = %f \n",second_order_z_term_i);

    //fprintf(screen,"third_order_x_term_i = %f \n",third_order_x_term_i);
    //fprintf(screen,"third_order_y_term_i = %f \n",third_order_y_term_i);
    //fprintf(screen,"third_order_z_term_i = %f \n",third_order_z_term_i);
  
    f[i][0] += gamma_factor*f_prefactor*(first_order_x_term[i] + second_order_x_term_i + third_order_x_term_i);
    f[i][1] += gamma_factor*f_prefactor*(first_order_y_term[i] + second_order_y_term_i + third_order_y_term_i);
    f[i][2] += gamma_factor*f_prefactor*(first_order_z_term[i] + second_order_z_term_i + third_order_z_term_i);

    //Bohm Potential calculation:

    bohm_pot = -gamma_factor*(f_prefactor/2)*((2*psi[i])/sum_packets[i] - (gamma_i)/sum_packets_i_2);
    if (eflag_global) eng_vdwl += bohm_pot;
    
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairBohmGlobalIons::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut,n+1,n+1,"pair:cut");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairBohmGlobalIons::settings(int narg, char **arg)
{
  if (narg != 5) error->all(FLERR,"Incorrect number of arguments!");

  cut_global = force->numeric(FLERR,arg[0]);
  sigma_global = force->numeric(FLERR,arg[1]);
  ion_sigma_global = force->numeric(FLERR,arg[2]);
  ion_core_electrons = force->numeric(FLERR,arg[3]);
  gamma_factor = force->numeric(FLERR,arg[4]);

  //fprintf(screen,"gamma_factor = %f \n",gamma_factor);
  //fprintf(screen,"ion_sigma_global = %f \n",ion_sigma_global);
  //fprintf(screen,"sigma_global = %f \n",sigma_global);

  hplanck  = force->hplanck;
  ////fprintf(screen,"hplanck = %f \n",hplanck);

  hbar = hplanck/(2*M_PI);

  e_sigma = sigma_global;
  e_sigma_4 = pow(e_sigma,4);
  e_sigma_3 = pow(e_sigma,3);
  e_sigma_2 = pow(e_sigma,2);

  ion_sigma = ion_sigma_global;
  ion_sigma_4 = pow(ion_sigma,4);
  ion_sigma_3 = pow(ion_sigma,3);
  ion_sigma_2 = pow(ion_sigma,2);

  e_packet_pre = 1/(e_sigma_3);
  ion_packet_pre = ion_core_electrons/(ion_sigma_3);

  //fprintf(screen,"e_packet_pre = %f \n",e_packet_pre);
  //fprintf(screen,"ion_packet_pre = %f \n",ion_packet_pre);
  
  ion_packet_zeta_sigma2 = ion_packet_pre/ion_sigma_2;
  e_packet_zeta_sigma2 = e_packet_pre/e_sigma_2;
  //fprintf(screen,"e_packet_zeta_sigma2 = %f \n",e_packet_zeta_sigma2);
  //fprintf(screen,"ion_packet_zeta_sigma2 = %f \n",ion_packet_zeta_sigma2);

  f_prefactor = (pow(hbar,2)/(4*e_mass));


  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairBohmGlobalIons::coeff(int narg, char **arg)
{
  if (narg < 2 || narg > 4)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

// Specify a different cutoff in pair_coeff step:

  double cut_one = cut_global;
  if (narg == 3) cut_one = force->numeric(FLERR,arg[2]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

// Specify a different sigma in pair_coeff step:

  double sigma_one = sigma_global;
  if (narg == 4) sigma_one = force->numeric(FLERR,arg[3]);

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}


/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairBohmGlobalIons::init_style()
{
  neighbor->request(this,instance_me);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairBohmGlobalIons::init_one(int i, int j)
{
  if (setflag[i][j] == 0)
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);

  return cut[i][j];
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairBohmGlobalIons::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) fwrite(&cut[i][j],sizeof(double),1,fp);
    }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairBohmGlobalIons::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR,&setflag[i][j],sizeof(int),1,fp,NULL,error);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) utils::sfread(FLERR,&cut[i][j],sizeof(double),1,fp,NULL,error);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairBohmGlobalIons::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairBohmGlobalIons::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR,&cut_global,sizeof(double),1,fp,NULL,error);
    utils::sfread(FLERR,&offset_flag,sizeof(int),1,fp,NULL,error);
    utils::sfread(FLERR,&mix_flag,sizeof(int),1,fp,NULL,error);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ---------------------------------------------------------------------- */

int PairBohmGlobalIons::pack_forward_comm(int n, int *list, double *buf,
                               int /*pbc_flag*/, int * /*pbc*/)
{
  int i,j,m;

  m = 0;
  if (commflag == 0){
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = sum_packets[j];
      buf[m++] = first_order_x_term[j];
      buf[m++] = first_order_y_term[j];
      buf[m++] = first_order_z_term[j];
      buf[m++] = theta_x[j];
      buf[m++] = theta_y[j];
      buf[m++] = theta_z[j];
      buf[m++] = psi[j];
      buf[m++] = rho_x[j];
      buf[m++] = rho_y[j];
      buf[m++] = rho_z[j];
      buf[m++] = rho_x1[j];
      buf[m++] = rho_x2[j];
      buf[m++] = rho_x3[j];
      buf[m++] = rho_y1[j];
      buf[m++] = rho_y3[j];
      buf[m++] = rho_z1[j];
    }
    return m;
  }
}

/* ---------------------------------------------------------------------- */

void PairBohmGlobalIons::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  if (commflag == 0){
    for (i = first; i < last; i++){

      sum_packets[i] = buf[m++];
      first_order_x_term[i] = buf[m++];
      first_order_y_term[i] = buf[m++];
      first_order_z_term[i] = buf[m++];

      theta_x[i] = buf[m++];
      theta_y[i] = buf[m++];
      theta_z[i] = buf[m++];

      psi[i] = buf[m++];

      rho_x[i] = buf[m++];
      rho_y[i] = buf[m++];
      rho_z[i] = buf[m++];

      rho_x1[i] = buf[m++];
      rho_x2[i] = buf[m++];
      rho_x3[i] = buf[m++];

      rho_y1[i] = buf[m++];
      rho_y3[i] = buf[m++];

      rho_z1[i] = buf[m++];
    }
  }
}

/* ---------------------------------------------------------------------- */

int PairBohmGlobalIons::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  if (commflag == 0){
    for (i = first; i < last; i++){
      buf[m++] = sum_packets[i];
      buf[m++] = first_order_x_term[i];
      buf[m++] = first_order_y_term[i];
      buf[m++] = first_order_z_term[i];

      buf[m++] = theta_x[i];
      buf[m++] = theta_y[i];
      buf[m++] = theta_z[i];

      buf[m++] = psi[i];

      buf[m++] = rho_x[i];
      buf[m++] = rho_y[i];
      buf[m++] = rho_z[i];

      buf[m++] = rho_x1[i];
      buf[m++] = rho_x2[i];
      buf[m++] = rho_x3[i];

      buf[m++] = rho_y1[i];
      buf[m++] = rho_y3[i];

      buf[m++] = rho_z1[i];
    }
    return m;
  }
}

/* ---------------------------------------------------------------------- */

void PairBohmGlobalIons::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  if (commflag == 0){
    for (i = 0; i < n; i++) {
      j = list[i];
      sum_packets[j] += buf[m++];
      first_order_x_term[j] += buf[m++];
      first_order_y_term[j] += buf[m++];
      first_order_z_term[j] += buf[m++];

      theta_x[j] += buf[m++];
      theta_y[j] += buf[m++];
      theta_z[j] += buf[m++];

      psi[j] += buf[m++];

      rho_x[j] += buf[m++];
      rho_y[j] += buf[m++];
      rho_z[j] += buf[m++];

      rho_x1[j] += buf[m++];
      rho_x2[j] += buf[m++];
      rho_x3[j] += buf[m++];

      rho_y1[j] += buf[m++];
      rho_y3[j] += buf[m++];

      rho_z1[j] += buf[m++];
    }
  }
}

/* ---------------------------------------------------------------------- */

void *PairBohmGlobalIons::extract(const char *str, int &dim)
{
  dim = 2;
  return NULL;
}
