/* -*- c++ -*- ----------------------------------------------------------
Bespoke pair_style to compute the BohmConserve force felt by electrons, according
to a Gaussian wave packet scheme.

Thomas Campbell
------------------------------------------------------------------------- */
#include "pair_bohm_conserve.h"
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

PairBohmConserve::PairBohmConserve(LAMMPS *lmp) : Pair(lmp) {
  nmax = 0;
  //fprintf(screen,"Using very very extra fast comm in pair_bohm_conserve...\n");
  manybody_flag = 1;

  sum_packets = NULL;
    
  first_order_x_term = NULL;
  first_order_y_term = NULL;
  first_order_z_term = NULL;
  
  second_order_x_term = NULL;
  second_order_y_term = NULL;
  second_order_z_term = NULL;

  third_order_x_term = NULL;
  third_order_y_term = NULL;
  third_order_z_term = NULL;
  
  theta_x = NULL;
  theta_y = NULL;
  theta_z = NULL;

  phi = NULL;

  eta = NULL;

  rho_x = NULL;
  rho_y = NULL;
  rho_z = NULL;
  
  rho_x1 = NULL;
  rho_x2 = NULL;
  rho_x3 = NULL;
  
  rho_y1 = NULL;
  rho_y3 = NULL;
  
  rho_z1 = NULL;

  first_order_sder = NULL;
  second_order_sder = NULL;
  third_order_sder = NULL;

  comm_forward = 27;
  comm_reverse = 27;
}

/* ---------------------------------------------------------------------- */

PairBohmConserve::~PairBohmConserve()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);

    memory->destroy(sum_packets);

    memory->destroy(first_order_x_term);
    memory->destroy(first_order_y_term);
    memory->destroy(first_order_z_term);
    
    memory->destroy(second_order_x_term);
    memory->destroy(second_order_y_term);
    memory->destroy(second_order_z_term);

    memory->destroy(third_order_x_term);
    memory->destroy(third_order_y_term);
    memory->destroy(third_order_z_term);
    
    memory->destroy(theta_x);
    memory->destroy(theta_y);
    memory->destroy(theta_z);

    memory->destroy(phi);

    memory->destroy(eta);

    memory->destroy(rho_x);
    memory->destroy(rho_y);
    memory->destroy(rho_z);
    
    memory->destroy(rho_x1);
    memory->destroy(rho_x2);
    memory->destroy(rho_x3);
    
    memory->destroy(rho_y1);
    memory->destroy(rho_y3);
    
    memory->destroy(rho_z1);

    memory->destroy(first_order_sder);
    memory->destroy(second_order_sder);
    memory->destroy(third_order_sder);

  }
}

/* ---------------------------------------------------------------------- */

void PairBohmConserve::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz;
  double xtmp2,ytmp2,ztmp2,delx2,dely2,delz2;
  double delx_2,dely_2,delz_2;
  double sigma_i,sigma_j;
  double sigma_i_7,sigma_i_5,sigma_i_2,sigma_i_3;
  double sigma_j_7,sigma_j_5,sigma_j_2;
  double sigma_i2,sigma_j2;
  double sigma_i2_6,sigma_i2_2,sigma_i2_5,sigma_i2_4;
  double sigma_j2_6,sigma_j2_2,sigma_j2_5,sigma_j2_4;
  double omega_ij;
  double omega_ji;
  double omega_ij2;
  double omega_ji2;
  double sum_packets_i_2;
  double sum_packets_i_3;
  double rsq;
  double rsq2;
  double reduced_theta_x,reduced_theta_y,reduced_theta_z;
  double reduced_theta_x_i,reduced_theta_y_i,reduced_theta_z_i;
  double beta_x,beta_y,beta_z;
  double beta_x_i,beta_y_i,beta_z_i;
  double first_self,second_self,third_self;
  double psi;
  double psi_i;
  double bohm_pot;
  double f_prefactor;
  double gamma_i;

  int *ilist,*jlist,*numneigh,**firstneigh;

  ev_init(eflag,vflag);

  if (atom->nmax > nmax) {
    memory->destroy(sum_packets);

    memory->destroy(first_order_x_term);
    memory->destroy(first_order_y_term);
    memory->destroy(first_order_z_term);
    
    memory->destroy(second_order_x_term);
    memory->destroy(second_order_y_term);
    memory->destroy(second_order_z_term);

    memory->destroy(third_order_x_term);
    memory->destroy(third_order_y_term);
    memory->destroy(third_order_z_term);
    
    memory->destroy(theta_x);
    memory->destroy(theta_y);
    memory->destroy(theta_z);

    memory->destroy(phi);

    memory->destroy(eta);

    memory->destroy(rho_x);
    memory->destroy(rho_y);
    memory->destroy(rho_z);
    
    memory->destroy(rho_x1);
    memory->destroy(rho_x2);
    memory->destroy(rho_x3);
    
    memory->destroy(rho_y1);
    memory->destroy(rho_y3);
    
    memory->destroy(rho_z1);

    memory->destroy(first_order_sder);
    memory->destroy(second_order_sder);
    memory->destroy(third_order_sder);

    nmax = atom->nmax;

    memory->create(sum_packets,nmax,"pair:sum_packets");

    memory->create(first_order_x_term,nmax,"pair:first_order_x_term");
    memory->create(first_order_y_term,nmax,"pair:first_order_y_term");
    memory->create(first_order_z_term,nmax,"pair:first_order_z_term");
    
    memory->create(second_order_x_term,nmax,"pair:second_order_x_term");
    memory->create(second_order_y_term,nmax,"pair:second_order_y_term");
    memory->create(second_order_z_term,nmax,"pair:second_order_z_term");

    memory->create(third_order_x_term,nmax,"pair:third_order_x_term");
    memory->create(third_order_y_term,nmax,"pair:third_order_y_term");
    memory->create(third_order_z_term,nmax,"pair:third_order_z_term");
    
    memory->create(theta_x,nmax,"pair:theta_x");
    memory->create(theta_y,nmax,"pair:theta_y");
    memory->create(theta_z,nmax,"pair:theta_z:");

    memory->create(phi,nmax,"pair:phi");

    memory->create(eta,nmax,"pair:eta");

    memory->create(rho_x,nmax,"pair:rho_x");
    memory->create(rho_y,nmax,"pair:rho_y");
    memory->create(rho_z,nmax,"pair:rho_z");
    
    memory->create(rho_x1,nmax,"pair:rho_x1");
    memory->create(rho_x2,nmax,"pair:rho_x2");
    memory->create(rho_x3,nmax,"pair:rho_x3");
    
    memory->create(rho_y1,nmax,"pair:rho_y1");
    memory->create(rho_y3,nmax,"pair:rho_y3");
    
    memory->create(rho_z1,nmax,"pair:rho_z1");

    memory->create(first_order_sder,nmax,"pair:first_order_sder");
    memory->create(second_order_sder,nmax,"pair:second_order_sder");
    memory->create(third_order_sder,nmax,"pair:third_order_sder");
  }

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;

  int *tag =atom->tag;

  double *mass = atom->mass;

  double *beta_correction = atom->beta_correction;

  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int newton_pair = force->newton_pair;

  double cutsquared = 0.;

  //hard code in electron mass for now

  double e_mass = 0.0005485799;

  hplanck  = force->hplanck;
  //fprintf(screen,"hplanck = %f \n",hplanck);

  hbar = hplanck/(2*M_PI);

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
      
      second_order_x_term[i] = 0.0;
      second_order_y_term[i] = 0.0;
      second_order_z_term[i] = 0.0;

      third_order_x_term[i] = 0.0;
      third_order_y_term[i] = 0.0;
      third_order_z_term[i] = 0.0;
      
      theta_x[i] = 0.0;
      theta_y[i] = 0.0;
      theta_z[i] = 0.0;

      phi[i] = 0.0;

      eta[i] = 0.0;
      rho_x[i] = 0.0;
      rho_y[i] = 0.0;
      rho_z[i] = 0.0;
      
      rho_x1[i] = 0.0;
      rho_x2[i] = 0.0;
      rho_x3[i] = 0.0;
      
      rho_y1[i] = 0.0;
      rho_y3[i] = 0.0;
      
      rho_z1[i] = 0.0;

      first_order_sder[i] = 0.0;
      second_order_sder[i] = 0.0;
      third_order_sder[i] = 0.0;
    }
  } 
  else{
    for (i = 0; i < nlocal; i++){
      sum_packets[i] = 0.0;
        
      first_order_x_term[i] = 0.0;
      first_order_y_term[i] = 0.0;
      first_order_z_term[i] = 0.0;
      
      second_order_x_term[i] = 0.0;
      second_order_y_term[i] = 0.0;
      second_order_z_term[i] = 0.0;

      third_order_x_term[i] = 0.0;
      third_order_y_term[i] = 0.0;
      third_order_z_term[i] = 0.0;
      
      theta_x[i] = 0.0;
      theta_y[i] = 0.0;
      theta_z[i] = 0.0;

      phi[i] = 0.0;

      eta[i] = 0.0;
      rho_x[i] = 0.0;
      rho_y[i] = 0.0;
      rho_z[i] = 0.0;
      
      rho_x1[i] = 0.0;
      rho_x2[i] = 0.0;
      rho_x3[i] = 0.0;
      
      rho_y1[i] = 0.0;
      rho_y3[i] = 0.0;
      
      rho_z1[i] = 0.0;

      second_order_sder[i] = 0.0;
      third_order_sder[i] = 0.0;
      first_order_sder[i] = 0.0;
    }
  }

  // loop over my atoms

  for (ii = 0; ii < inum; ii++) {

    i = ilist[ii];

    sigma_i = hplanck/(e_mass*sqrt(pow(v[i][0],2)+pow(v[i][1],2)+pow(v[i][2],2)));
    sigma_i_7 = pow(sigma_i,7);
    sigma_i_5 = pow(sigma_i,5);
    sigma_i_2 = pow(sigma_i,2);
    sigma_i_3 = pow(sigma_i,3);

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    itype = type[i];
    jlist = firstneigh[i];

    //fprintf(screen,"type[i] = %i \n",type[i]);

    jnum = numneigh[i];

     /* Self interaction terms */
    sum_packets[i] += 1/sigma_i_3;
    eta[i] += 3/sigma_i_5;
    phi[i] += -3/sigma_i_5;
    rho_x1[i] += -1/sigma_i_5;
    rho_y1[i] += -1/sigma_i_5;
    rho_z1[i] += -1/sigma_i_5;
    

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];

      rsq = delx*delx + dely*dely + delz*delz;

      jtype = type[j];

      //fprintf(screen,"type[j] = %i \n",type[j]);

      if (rsq < cutsquared) {

        delx_2 = pow(delx,2);
        dely_2 = pow(dely,2);
        delz_2 = pow(delz,2);

        sigma_j = hplanck/(e_mass*sqrt(pow(v[j][0],2)+pow(v[j][1],2)+pow(v[j][2],2)));
        sigma_j_7 = pow(sigma_j,7);
        sigma_j_5 = pow(sigma_j,5);
        sigma_j_2 = pow(sigma_j,2);

        omega_ij = exp(-(rsq)/(2*sigma_j_2));

        sum_packets[i] += omega_ij/(pow(sigma_j,3));
        first_order_x_term[i] += ((delx)/sigma_j_7)*(omega_ij)*(5 - (rsq)/sigma_j_2);
        first_order_y_term[i] += ((dely)/sigma_j_7)*(omega_ij)*(5 - (rsq)/sigma_j_2);
        first_order_z_term[i] += ((delz)/sigma_j_7)*(omega_ij)*(5 - (rsq)/sigma_j_2);
        theta_x[i] += ((delx)/sigma_j_5)*(omega_ij);
        theta_y[i] += ((dely)/sigma_j_5)*(omega_ij);
        theta_z[i] += ((delz)/sigma_j_5)*(omega_ij);
        phi[i] += ((rsq)/sigma_j_7 - 3/sigma_j_5)*omega_ij;
        eta[i] += (omega_ij)*(3/(sigma_j_5) - (rsq)/sigma_j_7);
        rho_x1[i] += (-1/sigma_j_5 + (delx_2)/sigma_j_7)*(omega_ij);
        rho_x2[i] += (((delx)*(dely))/sigma_j_7)*(omega_ij);
        rho_x3[i] += (((delx)*(delz))/sigma_j_7)*(omega_ij);
        rho_y1[i] += (-1/sigma_j_5 + (dely_2)/sigma_j_7)*(omega_ij);
        rho_y3[i] += (((dely)*(delz))/sigma_j_7)*(omega_ij);
        rho_z1[i] += (-1/sigma_j_5 + (delz_2)/sigma_j_7)*(omega_ij);

        if (newton_pair || j < nlocal) {

          omega_ji = exp(-(rsq)/(2*sigma_i_2));

          sum_packets[j] += omega_ji/sigma_i_3;
          first_order_x_term[j] -= ((delx)/sigma_i_7)*(omega_ji)*(5 - (rsq)/sigma_i_2);
          first_order_y_term[j] -= ((dely)/sigma_i_7)*(omega_ji)*(5 - (rsq)/sigma_i_2);
          first_order_z_term[j] -= ((delz)/sigma_i_7)*(omega_ji)*(5 - (rsq)/sigma_i_2);
          theta_x[j] -= ((delx)/sigma_i_5)*(omega_ji);
          theta_y[j] -= ((dely)/sigma_i_5)*(omega_ji);
          theta_z[j] -= ((delz)/sigma_i_5)*(omega_ji);
          phi[j] += ((rsq)/sigma_i_7 - 3/sigma_i_5)*omega_ji;
          eta[j] += (omega_ji)*(3/(sigma_i_5) - (rsq)/sigma_i_7);
          rho_x1[j] += (-1/sigma_i_5 + (delx_2)/sigma_i_7)*(omega_ji);
          rho_x2[j] += (((delx)*(dely))/sigma_i_7)*(omega_ji);
          rho_x3[j] += (((delx)*(delz))/sigma_i_7)*(omega_ji);
          rho_y1[j] += (-1/sigma_i_5 + (dely_2)/sigma_i_7)*(omega_ji);
          rho_y3[j] += (((dely)*(delz))/sigma_i_7)*(omega_ji);
          rho_z1[j] += (-1/sigma_i_5 + (delz_2)/sigma_i_7)*(omega_ji);

        }
      }
    }
  }
  comm_forward = 24;
  comm_reverse = 24;

  commflag = 0;

  if (newton_pair) comm->reverse_comm_pair(this);

  comm->forward_comm_pair(this);

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    
    sum_packets_i_2 = pow(sum_packets[i],2);
    sum_packets_i_3 = pow(sum_packets[i],3);

    first_order_x_term[i] *= (1/sum_packets[i]);
    first_order_y_term[i] *= (1/sum_packets[i]);
    first_order_z_term[i] *= (1/sum_packets[i]);

    rho_x[i] = theta_x[i]*rho_x1[i] + theta_y[i]*rho_x2[i] + theta_z[i]*rho_x3[i];
    
    rho_y[i] = theta_x[i]*rho_x2[i] + theta_y[i]*rho_y1[i] + theta_z[i]*rho_y3[i];

    rho_z[i] = theta_x[i]*rho_x3[i] + theta_y[i]*rho_y3[i] + theta_z[i]*rho_z1[i];

    second_order_x_term[i] = phi[i]*theta_x[i] + rho_x[i];
    second_order_y_term[i] = phi[i]*theta_y[i] + rho_y[i];
    second_order_z_term[i] = phi[i]*theta_z[i] + rho_z[i];
    
    second_order_x_term[i] *= (1/(sum_packets_i_2));
    second_order_y_term[i] *= (1/(sum_packets_i_2));
    second_order_z_term[i] *= (1/(sum_packets_i_2));

    gamma_i = pow(theta_x[i],2) + pow(theta_y[i],2) + pow(theta_z[i],2);
    
    third_order_x_term[i] = -theta_x[i]*gamma_i;
    third_order_y_term[i] = -theta_y[i]*gamma_i;
    third_order_z_term[i] = -theta_z[i]*gamma_i;
    
    third_order_x_term[i] *= (1/(sum_packets_i_3));
    third_order_y_term[i] *= (1/(sum_packets_i_3));
    third_order_z_term[i] *= (1/(sum_packets_i_3));

    f_prefactor = (pow(hbar,2)/(4*e_mass));
    
    f[i][0] += f_prefactor*(first_order_x_term[i] + second_order_x_term[i] + third_order_x_term[i]);
    f[i][1] += f_prefactor*(first_order_y_term[i] + second_order_y_term[i] + third_order_y_term[i]);
    f[i][2] += f_prefactor*(first_order_z_term[i] + second_order_z_term[i] + third_order_z_term[i]);

    //Bohm Potential calculation:

    bohm_pot = -(f_prefactor/2)*((-2*eta[i])/sum_packets[i] - (gamma_i)/sum_packets_i_2);
    if (eflag_global) eng_vdwl += bohm_pot;
    
    sigma_i2 = hplanck/(e_mass*sqrt(pow(v[i][0],2)+pow(v[i][1],2)+pow(v[i][2],2)));

    xtmp2 = x[i][0];
    ytmp2 = x[i][1];
    ztmp2 = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    /* Add self interaction terms */

    sigma_i2_2 = pow(sigma_i2,2);
    sigma_i2_6 = pow(sigma_i2,6);
    sigma_i2_5 = pow(sigma_i2,5);
    sigma_i2_4 = pow(sigma_i2,4);


    first_order_sder[i] += (1/sum_packets[i])*(30/sigma_i2_6);
    first_self = (1/sum_packets[i])*(30/sigma_i2_6);

    second_order_sder[i] += (1/sum_packets_i_2)*(2*eta[i]*(-3/sigma_i2_4));
    second_self = (1/sum_packets_i_2)*(2*eta[i]*(-3/sigma_i2_4));

    third_order_sder[i] += (1/sum_packets_i_3)*(2*gamma_i*(-3/sigma_i2_4));
    third_self = (1/sum_packets_i_3)*(2*gamma_i*(-3/sigma_i2_4));

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      jtype = type[j];

      /* Careful with indices. Note that these 'del' terms are j - i. */

      delx2 = - xtmp2 + x[j][0];
      dely2 = - ytmp2 + x[j][1];
      delz2 = - ztmp2 + x[j][2];

      rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;

      if (rsq2 < cutsquared) {

        omega_ji2 = exp(-(rsq2)/(2*sigma_i2_2));

        first_order_sder[i] += (1/sum_packets[j])*(30/(sigma_i2_6) - (20*rsq2)/(pow(sigma_i2,8)) + (2*pow(rsq2,2))/(pow(sigma_i2,10)))*omega_ji2;
        reduced_theta_x_i = theta_x[j] - (delx2*omega_ji2/sigma_i2_5);
        reduced_theta_y_i = theta_y[j] - (dely2*omega_ji2/sigma_i2_5);
        reduced_theta_z_i = theta_z[j] - (delz2*omega_ji2/sigma_i2_5);
        beta_x_i = (delx2*omega_ji2/sigma_i2_6)*(5 - rsq2/sigma_i2_2);
        beta_y_i = (dely2*omega_ji2/sigma_i2_6)*(5 - rsq2/sigma_i2_2);
        beta_z_i = (delz2*omega_ji2/sigma_i2_6)*(5 - rsq2/sigma_i2_2);
        psi_i = 2*(beta_x_i*reduced_theta_x_i + beta_y_i*reduced_theta_y_i + beta_z_i*reduced_theta_z_i);
        second_order_sder[i] += (1/pow(sum_packets[j],2))*(2*eta[j]*(rsq2/sigma_i2_2 - 3)*(omega_ji2/sigma_i2_4) + (10 - (2*rsq2)/sigma_i2_2)*(rsq2*pow(omega_ji2,2)/pow(sigma_i2,11)) + psi_i);
        third_order_sder[i] += (1/pow(sum_packets[j],3))*((pow(theta_x[j],2) + pow(theta_y[j],2) + pow(theta_z[j],2))*2*(rsq2/sigma_i2_2 - 3)*(omega_ji2/sigma_i2_4));

        if (newton_pair || j < nlocal) {

          sigma_j2 = hplanck/(e_mass*sqrt(pow(v[j][0],2)+pow(v[j][1],2)+pow(v[j][2],2)));

          sigma_j2_2 = pow(sigma_j2,2);
          sigma_j2_6 = pow(sigma_j2,6);
          sigma_j2_5 = pow(sigma_j2,5);
          sigma_j2_4 = pow(sigma_j2,4);

          omega_ij2 = exp(-(rsq2)/(2*sigma_j2_2));

          first_order_sder[j] += (1/sum_packets[i])*(30/(sigma_j2_6) - (20*rsq2)/(pow(sigma_j2,8)) + (2*pow(rsq2,2))/(pow(sigma_j2,10)))*omega_ij2;
          reduced_theta_x = theta_x[i] + (delx2*omega_ij2/sigma_j2_5);
          reduced_theta_y = theta_y[i] + (dely2*omega_ij2/sigma_j2_5);
          reduced_theta_z = theta_z[i] + (delz2*omega_ij2/sigma_j2_5);
          beta_x = - (delx2*omega_ij2/sigma_j2_6)*(5 - rsq2/sigma_j2_2);
          beta_y = - (dely2*omega_ij2/sigma_j2_6)*(5 - rsq2/sigma_j2_2);
          beta_z = - (delz2*omega_ij2/sigma_j2_6)*(5 - rsq2/sigma_j2_2);
          psi = 2*(beta_x*reduced_theta_x + beta_y*reduced_theta_y + beta_z*reduced_theta_z);
          second_order_sder[j] += (1/sum_packets_i_2)*(2*eta[i]*(rsq2/sigma_j2_2 - 3)*(omega_ij2/sigma_j2_4) + (10 - (2*rsq2)/sigma_j2_2)*(rsq2*pow(omega_ij2,2)/pow(sigma_j2,11)) + psi);
          third_order_sder[j] += (1/sum_packets_i_3)*((gamma_i)*2*(rsq2/sigma_j2_2 - 3)*(omega_ij2/sigma_j2_4));
        }
      }
    }
  }
  commflag = 1;

  comm_reverse = 3;

  if (newton_pair) comm->reverse_comm_pair(this);

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    beta_correction[i] = (hplanck/(pow(e_mass,2)*pow((pow(v[i][0],2)+pow(v[i][1],2)+pow(v[i][2],2)),1.5)))*(-(f_prefactor/2)*(first_order_sder[i] + second_order_sder[i] + third_order_sder[i]));
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairBohmConserve::allocate()
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

void PairBohmConserve::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal pair_style command");

  cut_global = force->numeric(FLERR,arg[0]);

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

void PairBohmConserve::coeff(int narg, char **arg)
{
  if (narg < 2 || narg > 3)
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

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}


/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairBohmConserve::init_style()
{
  if (!atom->beta_flag)
    error->all(FLERR,"Pair style bohm_conserve requires atom attribute beta");
  neighbor->request(this,instance_me);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairBohmConserve::init_one(int i, int j)
{
  if (setflag[i][j] == 0)
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);

  return cut[i][j];
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairBohmConserve::write_restart(FILE *fp)
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

void PairBohmConserve::read_restart(FILE *fp)
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

void PairBohmConserve::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairBohmConserve::read_restart_settings(FILE *fp)
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

int PairBohmConserve::pack_forward_comm(int n, int *list, double *buf,
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
      buf[m++] = second_order_x_term[j];
      buf[m++] = second_order_y_term[j];
      buf[m++] = second_order_z_term[j];
      buf[m++] = third_order_x_term[j];
      buf[m++] = third_order_y_term[j];
      buf[m++] = third_order_z_term[j];
      buf[m++] = theta_x[j];
      buf[m++] = theta_y[j];
      buf[m++] = theta_z[j];
      buf[m++] = phi[j];
      buf[m++] = eta[j];
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
  if (commflag == 1){
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = first_order_sder[j];
      buf[m++] = second_order_sder[j];
      buf[m++] = third_order_sder[j];
    }
    return m;
  }
}

/* ---------------------------------------------------------------------- */

void PairBohmConserve::unpack_forward_comm(int n, int first, double *buf)
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

      second_order_x_term[i] = buf[m++];
      second_order_y_term[i] = buf[m++];
      second_order_z_term[i] = buf[m++];

      third_order_x_term[i] = buf[m++];
      third_order_y_term[i] = buf[m++];
      third_order_z_term[i] = buf[m++];

      theta_x[i] = buf[m++];
      theta_y[i] = buf[m++];
      theta_z[i] = buf[m++];

      phi[i] = buf[m++];

      eta[i] = buf[m++];

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
  if (commflag == 1){
    for (i = first; i < last; i++){
      first_order_sder[i] = buf[m++];
      second_order_sder[i] = buf[m++];
      third_order_sder[i] = buf[m++];
    }
  }
}

/* ---------------------------------------------------------------------- */

int PairBohmConserve::pack_reverse_comm(int n, int first, double *buf)
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

      buf[m++] = second_order_x_term[i];
      buf[m++] = second_order_y_term[i];
      buf[m++] = second_order_z_term[i];

      buf[m++] = third_order_x_term[i];
      buf[m++] = third_order_y_term[i];
      buf[m++] = third_order_z_term[i];

      buf[m++] = theta_x[i];
      buf[m++] = theta_y[i];
      buf[m++] = theta_z[i];

      buf[m++] = phi[i];

      buf[m++] = eta[i];

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
  if (commflag == 1){
    for (i = first; i < last; i++){
      buf[m++] = first_order_sder[i];
      buf[m++] = second_order_sder[i];
      buf[m++] = third_order_sder[i];
    }
    return m;
  }
}

/* ---------------------------------------------------------------------- */

void PairBohmConserve::unpack_reverse_comm(int n, int *list, double *buf)
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

      second_order_x_term[j] += buf[m++];
      second_order_y_term[j] += buf[m++];
      second_order_z_term[j] += buf[m++];

      third_order_x_term[j] += buf[m++];
      third_order_y_term[j] += buf[m++];
      third_order_z_term[j] += buf[m++];

      theta_x[j] += buf[m++];
      theta_y[j] += buf[m++];
      theta_z[j] += buf[m++];

      phi[j] += buf[m++];

      eta[j] += buf[m++];

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
  if (commflag == 1){
    for (i = 0; i < n; i++) {
      j = list[i];
      first_order_sder[j] += buf[m++];
      second_order_sder[j] += buf[m++];
      third_order_sder[j] += buf[m++];
    }
  }
}

/* ---------------------------------------------------------------------- */

void *PairBohmConserve::extract(const char *str, int &dim)
{
  dim = 2;
  return NULL;
}
