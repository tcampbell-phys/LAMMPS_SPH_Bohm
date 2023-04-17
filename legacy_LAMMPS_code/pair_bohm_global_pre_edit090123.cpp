/* -*- c++ -*- ----------------------------------------------------------
Bespoke pair_style to compute the BohmGlobal force felt by electrons, according
to a Gaussian wave packet scheme.

Thomas Campbell
------------------------------------------------------------------------- */
#include "pair_bohm_global.h"
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

PairBohmGlobal::PairBohmGlobal(LAMMPS *lmp) : Pair(lmp) {
  nmax = 0;
  //fprintf(screen,"Using very very extra fast comm in pair_bohm_Global...\n");
  manybody_flag = 1;

  sum_packets = NULL;
    
  first_order_x_term = NULL;
  first_order_y_term = NULL;
  first_order_z_term = NULL;
  
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

  comm_forward = 18;
  comm_reverse = 18;
}

/* ---------------------------------------------------------------------- */

PairBohmGlobal::~PairBohmGlobal()
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
  }
}

/* ---------------------------------------------------------------------- */

void PairBohmGlobal::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz;
  double xtmp2,ytmp2,ztmp2,delx2,dely2,delz2;
  double delx_2,dely_2,delz_2;
  double sigma;
  double sigma_4,sigma_2;
  double omega_ij;
  double omega_ij2;
  double sum_packets_i_2;
  double sum_packets_i_3;
  double rsq;
  double rsq2;
  double reduced_theta_x,reduced_theta_y,reduced_theta_z;
  double reduced_theta_x_i,reduced_theta_y_i,reduced_theta_z_i;
  double bohm_pot;
  double f_prefactor;
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

    nmax = atom->nmax;

    memory->create(sum_packets,nmax,"pair:sum_packets");

    memory->create(first_order_x_term,nmax,"pair:first_order_x_term");
    memory->create(first_order_y_term,nmax,"pair:first_order_y_term");
    memory->create(first_order_z_term,nmax,"pair:first_order_z_term");
    
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
  }

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;

  if (atom->forceterms_flag == 1){
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
  }

  int *tag =atom->tag;

  double *mass = atom->mass;

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
    }
  }

  // loop over my atoms

  sigma = sigma_global;
  //fprintf(screen,"sigma_global = %f \n",sigma_global);
  sigma_4 = pow(sigma,4);
  sigma_2 = pow(sigma,2);

  /* Note conditions on type in below loops. These are applied to enable Ions (type 1) to experience the Bohm
  field but not to contribute. Have checked that the electron Bohm forces are not affected. */

  for (ii = 0; ii < inum; ii++) {

    i = ilist[ii];

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    itype = type[i];
    jlist = firstneigh[i];

    //fprintf(screen,"type[i] = %i \n",type[i]);

    jnum = numneigh[i];

     /* Self interaction terms */
    sum_packets[i] += 1;
    eta[i] += 3/sigma_2;
    phi[i] += -3/sigma_2;
    rho_x1[i] += -1/sigma_2;
    rho_y1[i] += -1/sigma_2;
    rho_z1[i] += -1/sigma_2;
    

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

        omega_ij = exp(-(rsq)/(2*sigma_2));

        sum_packets[i] += omega_ij;
        first_order_x_term[i] += ((delx)/sigma_4)*(omega_ij)*(5 - (rsq)/sigma_2);
        first_order_y_term[i] += ((dely)/sigma_4)*(omega_ij)*(5 - (rsq)/sigma_2);
        first_order_z_term[i] += ((delz)/sigma_4)*(omega_ij)*(5 - (rsq)/sigma_2);
        theta_x[i] += ((delx)/sigma_2)*(omega_ij);
        theta_y[i] += ((dely)/sigma_2)*(omega_ij);
        theta_z[i] += ((delz)/sigma_2)*(omega_ij);
        phi[i] += ((rsq)/sigma_4 - 3/sigma_2)*omega_ij;
        eta[i] += (omega_ij)*(3/(sigma_2) - (rsq)/sigma_4);
        rho_x1[i] += (-1/sigma_2 + (delx_2)/sigma_4)*(omega_ij);
        rho_x2[i] += (((delx)*(dely))/sigma_4)*(omega_ij);
        rho_x3[i] += (((delx)*(delz))/sigma_4)*(omega_ij);
        rho_y1[i] += (-1/sigma_2 + (dely_2)/sigma_4)*(omega_ij);
        rho_y3[i] += (((dely)*(delz))/sigma_4)*(omega_ij);
        rho_z1[i] += (-1/sigma_2 + (delz_2)/sigma_4)*(omega_ij);


        if (newton_pair || j < nlocal) {

          sum_packets[j] += omega_ij;
          first_order_x_term[j] -= ((delx)/sigma_4)*(omega_ij)*(5 - (rsq)/sigma_2);
          first_order_y_term[j] -= ((dely)/sigma_4)*(omega_ij)*(5 - (rsq)/sigma_2);
          first_order_z_term[j] -= ((delz)/sigma_4)*(omega_ij)*(5 - (rsq)/sigma_2);
          theta_x[j] -= ((delx)/sigma_2)*(omega_ij);
          theta_y[j] -= ((dely)/sigma_2)*(omega_ij);
          theta_z[j] -= ((delz)/sigma_2)*(omega_ij);
          phi[j] += ((rsq)/sigma_4 - 3/sigma_2)*omega_ij;
          eta[j] += (omega_ij)*(3/(sigma_2) - (rsq)/sigma_4);
          rho_x1[j] += (-1/sigma_2 + (delx_2)/sigma_4)*(omega_ij);
          rho_x2[j] += (((delx)*(dely))/sigma_4)*(omega_ij);
          rho_x3[j] += (((delx)*(delz))/sigma_4)*(omega_ij);
          rho_y1[j] += (-1/sigma_2 + (dely_2)/sigma_4)*(omega_ij);
          rho_y3[j] += (((dely)*(delz))/sigma_4)*(omega_ij);
          rho_z1[j] += (-1/sigma_2 + (delz_2)/sigma_4)*(omega_ij);

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

    //fprintf(screen,"sum_packets[i] = %20.15g \n",sum_packets[i]);

    if (atom->forceterms_flag == 1){
      sum_packets_d[i] = sum_packets[i];
      first_order_x_term_d[i] = first_order_x_term[i];
      first_order_y_term_d[i] = first_order_y_term[i];
      first_order_z_term_d[i] = first_order_z_term[i];
    }
    //fprintf(screen,"first_order_x_term[i] = %20.15g \n",first_order_x_term[i]);
    //fprintf(screen,"first_order_y_term[i] = %20.15g \n",first_order_y_term[i]);
    //fprintf(screen,"first_order_z_term[i] = %20.15g \n",first_order_z_term[i]);

    first_order_x_term[i] *= (1/sum_packets[i]);
    first_order_y_term[i] *= (1/sum_packets[i]);
    first_order_z_term[i] *= (1/sum_packets[i]);

    rho_x[i] = theta_x[i]*rho_x1[i] + theta_y[i]*rho_x2[i] + theta_z[i]*rho_x3[i];
    
    rho_y[i] = theta_x[i]*rho_x2[i] + theta_y[i]*rho_y1[i] + theta_z[i]*rho_y3[i];

    rho_z[i] = theta_x[i]*rho_x3[i] + theta_y[i]*rho_y3[i] + theta_z[i]*rho_z1[i];

    //fprintf(screen,"rho_x[i] = %20.15g \n",rho_x[i]);
    //fprintf(screen,"rho_y[i] = %20.15g \n",rho_y[i]);
    //fprintf(screen,"rho_z[i] = %20.15g \n",rho_z[i]);

    //fprintf(screen,"phi[i] = %20.15g \n",phi[i]);

    //fprintf(screen,"theta_x[i] = %20.15g \n",theta_x[i]);
    //fprintf(screen,"theta_y[i] = %20.15g \n",theta_y[i]);
    //fprintf(screen,"theta_z[i] = %20.15g \n",theta_z[i]);

    second_order_x_term_i = phi[i]*theta_x[i] + rho_x[i];
    second_order_y_term_i = phi[i]*theta_y[i] + rho_y[i];
    second_order_z_term_i = phi[i]*theta_z[i] + rho_z[i];

    if (atom->forceterms_flag == 1){
      phi_d[i] = phi[i];

      rho_x_d[i] = rho_x[i];
      rho_y_d[i] = rho_z[i];
      rho_z_d[i] = rho_x[i];

      theta_x_d[i] = theta_x[i];
      theta_y_d[i] = theta_z[i];
      theta_z_d[i] = theta_x[i];

      second_order_x_term_d[i] = second_order_x_term_i;
      second_order_y_term_d[i] = second_order_y_term_i;
      second_order_z_term_d[i] = second_order_z_term_i;
    }

    //fprintf(screen,"second_order_x_term_i = %20.15g \n",second_order_x_term_i);
    //fprintf(screen,"second_order_y_term_i = %20.15g \n",second_order_y_term_i);
    //fprintf(screen,"second_order_z_term_i = %20.15g \n",second_order_z_term_i);

    
    second_order_x_term_i *= (1/(sum_packets_i_2));
    second_order_y_term_i *= (1/(sum_packets_i_2));
    second_order_z_term_i *= (1/(sum_packets_i_2));

    gamma_i = pow(theta_x[i],2) + pow(theta_y[i],2) + pow(theta_z[i],2);
    
    third_order_x_term_i = -theta_x[i]*gamma_i;
    third_order_y_term_i = -theta_y[i]*gamma_i;
    third_order_z_term_i = -theta_z[i]*gamma_i;

    //fprintf(screen,"third_order_x_term_i = %20.15g \n",third_order_x_term_i);
    //fprintf(screen,"third_order_y_term_i = %20.15g \n",third_order_y_term_i);
    //fprintf(screen,"third_order_z_term_i = %20.15g \n",third_order_z_term_i);

    if (atom->forceterms_flag == 1){
      third_order_x_term_d[i] = third_order_x_term_i;
      third_order_y_term_d[i] = third_order_y_term_i;
      third_order_z_term_d[i] = third_order_z_term_i;
    }
    
    third_order_x_term_i *= (1/(sum_packets_i_3));
    third_order_y_term_i *= (1/(sum_packets_i_3));
    third_order_z_term_i *= (1/(sum_packets_i_3));

    f_prefactor = (pow(hbar,2)/(4*e_mass));
    
    f[i][0] += gamma_factor*f_prefactor*(first_order_x_term[i] + second_order_x_term_i + third_order_x_term_i);
    f[i][1] += gamma_factor*f_prefactor*(first_order_y_term[i] + second_order_y_term_i + third_order_y_term_i);
    f[i][2] += gamma_factor*f_prefactor*(first_order_z_term[i] + second_order_z_term_i + third_order_z_term_i);

    //Bohm Potential calculation:

    bohm_pot = -gamma_factor*(f_prefactor/2)*((-2*eta[i])/sum_packets[i] - (gamma_i)/sum_packets_i_2);
    if (eflag_global) eng_vdwl += bohm_pot;
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairBohmGlobal::allocate()
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

void PairBohmGlobal::settings(int narg, char **arg)
{
  if (narg != 3) error->all(FLERR,"Illegal pair_style command");

  cut_global = force->numeric(FLERR,arg[0]);
  sigma_global = force->numeric(FLERR,arg[1]);
  gamma_factor = force->numeric(FLERR,arg[2]);

  fprintf(screen,"gamma_factor = %f \n",gamma_factor);

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

void PairBohmGlobal::coeff(int narg, char **arg)
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

void PairBohmGlobal::init_style()
{
  neighbor->request(this,instance_me);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairBohmGlobal::init_one(int i, int j)
{
  if (setflag[i][j] == 0)
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);

  return cut[i][j];
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairBohmGlobal::write_restart(FILE *fp)
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

void PairBohmGlobal::read_restart(FILE *fp)
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

void PairBohmGlobal::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairBohmGlobal::read_restart_settings(FILE *fp)
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

int PairBohmGlobal::pack_forward_comm(int n, int *list, double *buf,
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
}

/* ---------------------------------------------------------------------- */

void PairBohmGlobal::unpack_forward_comm(int n, int first, double *buf)
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
}

/* ---------------------------------------------------------------------- */

int PairBohmGlobal::pack_reverse_comm(int n, int first, double *buf)
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
}

/* ---------------------------------------------------------------------- */

void PairBohmGlobal::unpack_reverse_comm(int n, int *list, double *buf)
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
}

/* ---------------------------------------------------------------------- */

void *PairBohmGlobal::extract(const char *str, int &dim)
{
  dim = 2;
  return NULL;
}
