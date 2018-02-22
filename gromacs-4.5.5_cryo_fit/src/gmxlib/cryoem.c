/* Cryo-EM fitting
   reference: Flexible Multi-scale Fitting of Atomic Structures into Low-resolution Electron Density Maps with Elastic Network Normal Mode Analysis
              Florence Tama, JMB, 2004
*/
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cryoem.h"
#include "smalloc.h"
#include "gmx_fatal.h"
#include "network.h"
#include <math.h>
static real em_consErf=0.0, em_consExp=0.0;
static real em_pow1=0.0, em_pow2=0.0;
static real PI=0.0;
static int vsize=0;

static void read_emmap(char *emfile, t_emmap *emmap)
{
//  printf("start read_emmap\n");

  FILE *fp;
  int stat,i,j,k;
  int vol=0;
  if((fp=fopen(emfile,"r")) == NULL) {
    gmx_fatal(FARGS, "Cannot open emmap file!");
  }

//  stat=fscanf(fp,"%f %f %f %f %d %d %d %f\n",&(emmap->dx),&(emmap->x0),&(emmap->y0),&(emmap->z0),&(emmap->nx),&(emmap->ny),&(emmap->nz),&(emmap->thresh));
  stat=fscanf(fp,"%f %f %f %f %d %d %d\n",&(emmap->dx),&(emmap->x0),&(emmap->y0),&(emmap->z0),&(emmap->nx),&(emmap->ny),&(emmap->nz));
//printf("Voxels with values less than %f will be set to %f\n",emmap->thresh, emmap->thresh);	
  if (stat != 7) {
    gmx_fatal(FARGS,"Error reading EM density map!");
  }
  emmap->dx=emmap->dx/10.0;
  emmap->dy=emmap->dx;
  emmap->dz=emmap->dx;
  
  emmap->x0=emmap->x0/10.0;
  emmap->y0=emmap->y0/10.0;
  emmap->z0=emmap->z0/10.0;
  
  vol = emmap->nx * emmap->ny * emmap->nz;
  snew(emmap->rho_em, vol);
/* TRIPLE */    
  for(k=0;k<emmap->nz;k++)
    for(j=0;j<emmap->ny;j++)
      for(i=0;i<emmap->nx;i++) {
        vol = i*emmap->ny*emmap->nz + j*emmap->nz + k;
        fscanf(fp,"%f ",&emmap->rho_em[vol]);

/* THRESHOLD */
	if(emmap->rho_em[vol] < emmap->thresh){
	emmap->rho_em[vol]=0;
	} 
	
      }
  fclose(fp);

/* done reading in the em map, now read in the variance map. */

}

/* Calculate the cutoff distance for density evaluation */
static void setCutoff(t_emmap *simmap)
{
  real r=0.0, y;
  real cons1;
  
  cons1=sqrt(6.0/PI);
  while ( TRUE ) {
    y = erf(sqrt(1.5*r)) - cons1*r*exp(-1.5*r*r);
    
    if (1-y < simmap->err)
      break;
    r += 0.01;
  }
  
  simmap->nGridCutX = (int)ceilf(r*simmap->sigma*simmap->rdx);
  simmap->nGridCutY = (int)ceilf(r*simmap->sigma*simmap->rdy);
  simmap->nGridCutZ = (int)ceilf(r*simmap->sigma*simmap->rdz);
  
}

static void collect_emmap_rho(real *rho, int size, t_commrec *cr)
{

  real rho_tot[size];

#ifdef GMX_MPI
  MPI_Allreduce(rho, rho_tot, size, MPI_FLOAT, MPI_SUM, cr->mpi_comm_mysim); /*collect partial rho(s)*/
  memcpy(rho, rho_tot, size*sizeof(float));
#endif
}

static void bound_check(t_mdatoms *md, rvec x[], real lx, real hx, real ly, real hy, real lz, real hz) 
{
  int i, start, end;
  start=md->start;
  end=md->start+md->homenr;
  
  for (i=start; i< end; i++) {
    if (x[i][0]<lx || x[i][0]>hx)
      gmx_fatal(FARGS, "x is beyond the map boundary!");
    if (x[i][1]<ly || x[i][1]>hy)
      gmx_fatal(FARGS, "y is beyond the map boundary!");
    if (x[i][2]<lz || x[i][2]>hz)
      gmx_fatal(FARGS, "y is beyond the map boundary!");
  }

}

static void emmap_zero(t_emmap *emmap)
{
  int i,j,k,n, index;

/* TRIPLE */
  for(i=0;i<emmap->nx;i++)
    for(j=0;j<emmap->ny;j++)
      for(k=0;k<emmap->nz;k++) {
        index = i* emmap->ny * emmap->nz + j *emmap->nz +k;
	emmap->rho_sim[index]=0.0;
      }
}
static void ig_check(t_emmap *emmap, gmx_groups_t *group, int mdfitgrp, t_commrec *cr, t_state *state)
{
  int n, start, end;
  int *igx, *igy, *igz;
  char info[512];
  start=0;
  end=state->natoms;
   
  igx=emmap->ig_x;
  igy=emmap->ig_y;
  igz=emmap->ig_z;
   
  for(n=start;n<end;n++) {
    if (ggrpnr(group,egcMDFIT,n) < mdfitgrp) {
      if (igx[n]-emmap->nGridCutX < 0 ) {
      igx[n] = emmap->nGridCutX;
      }
      if (igx[n]+emmap->nGridCutX >= emmap->nx) {
         igx[n] = emmap->nx - emmap->nGridCutX - 1;
      }
      if (igy[n]-emmap->nGridCutY < 0 ) {
        igy[n] = emmap->nGridCutY;
      }
      if ( igy[n]+emmap->nGridCutY >= emmap->ny) {
          igy[n] = emmap->ny - emmap->nGridCutY -1;
      }      
      if (igz[n]-emmap->nGridCutZ < 0 ){
          igz[n]=emmap->nGridCutZ;
      }
      if (igz[n]+emmap->nGridCutZ >= emmap->nz) {
         igz[n]=emmap->nz-emmap->nGridCutZ-1;
      }
    }
  }
}

      
static void ig_check1(t_emmap *emmap, t_mdatoms *md, int mdfitgrp, t_commrec *cr, rvec x[])
{
  int n, start, end;
  int *igx, *igy, *igz;
  char info[512];
  
  start=md->start;
  end=md->start+md->homenr;
   
  igx=emmap->ig_x;
  igy=emmap->ig_y;
  igz=emmap->ig_z;
   
  for(n=start;n<end;n++) {
    if (md->cMDFIT[n] < mdfitgrp) {
      if (igx[n]-emmap->nGridCutX < 0 || igx[n]+emmap->nGridCutX >= emmap->nx) {
        sprintf(info,"Atom %d on grid x:%d (%f) is out of map!",n, igx[n],x[n][0]);
        gmx_fatal(FARGS,info);
      }
      if (igy[n]-emmap->nGridCutY < 0 || igy[n]+emmap->nGridCutY >= emmap->ny) {
        sprintf(info,"Atom %d on grid y:%d (%f) is out of map!", n, igy[n],x[n][1]);
        gmx_fatal(FARGS,info); 
      }
      if (igz[n]-emmap->nGridCutZ < 0 || igz[n]+emmap->nGridCutZ >= emmap->nz) {
        sprintf(info,"Atom %d on grid z:%d (%f) is out of map!", n, igz[n],x[n][2]);
        gmx_fatal(FARGS,info);
      }
    }
  }
}

static real corrCoef(t_emmap *emmap)
{
  real cc=0.0, norm_exp=0.0, norm_sim=0.0, exp_sim=0.0;
  int i,j,k, index;
/* TRIPLE */
  for(i=0;i<emmap->nx;i++) 
    for(j=0;j<emmap->ny;j++)
      for(k=0;k<emmap->nz;k++) {
/* his may be the problem line*/
         index = i* emmap->ny * emmap->nz + j *emmap->nz +k;
         norm_sim += emmap->rho_sim[index]*emmap->rho_sim[index];
/* paul addition*/
	 exp_sim += emmap->rho_em[index]*emmap->rho_sim[index];
      }

  norm_exp=sqrt(emmap->rho2_em);
  norm_sim=sqrt(norm_sim);
  
  cc=exp_sim/(norm_exp*norm_sim);
  
  return cc;
}

void expand_emmap(t_emmap *emmap, t_state *state, t_commrec *cr)
{
//  printf("start expand_emmap\n");
  int nr;
  nr = state->natoms;
  srenew(emmap->ig_x,nr);
  srenew(emmap->ig_y,nr);
  srenew(emmap->ig_z,nr);
  srenew(emmap->derfa_x, nr*emmap->nx);
  srenew(emmap->derfa_y, nr*emmap->ny);
  srenew(emmap->derfa_z, nr*emmap->nz);
  srenew(emmap->dexpa_x, nr*emmap->nx);
  srenew(emmap->dexpa_y, nr*emmap->ny);
  srenew(emmap->dexpa_z, nr*emmap->nz);
  memset(emmap->derfa_x,0,nr*emmap->nx);
  memset(emmap->derfa_y,0,nr*emmap->ny);
  memset(emmap->derfa_z,0,nr*emmap->nz);
  memset(emmap->dexpa_x,0,nr*emmap->nx);
  memset(emmap->dexpa_y,0,nr*emmap->ny);
  memset(emmap->dexpa_z,0,nr*emmap->nz);
    
  memset(emmap->ig_x,0,nr);
  memset(emmap->ig_y,0,nr);
  memset(emmap->ig_z,0,nr);

  emmap->nr = nr;
}

/* initialization of *emmap and *simmap */       
void init_em_map(FILE *fplog, char *emfile, t_emmap *emmap, t_inputrec *ir, t_state *state,  t_commrec *cr)
{
   int n,i,j,k,index;
   double norm_exp=0.0,sigma;
   double ave_exp=0.0;
   int nr;
   nr=state->natoms;
   emmap->sigma = ir->emsigma;
   emmap->err = ir->emcutoff;
   emmap->steps = ir->emsteps;
   emmap->weight = ir->emweight;
   emmap->emwritefrequency = ir->emwritefrequency;
   emmap->thresh = ir->emthresh;
  /* constants used in rho and rho' evaluations */
   sigma=emmap->sigma;
   em_consErf=sqrt(3/(2*sigma*sigma));
   em_consExp=-3.0/(2*sigma*sigma);
   PI=acos(-1);
   em_pow1=1.0/sqrt(PI*sigma*sigma/6.0);
   emmap->rho_em = NULL;
   emmap->rho_sim = NULL;
   emmap->derfa_x = NULL;
   emmap->derfa_y = NULL;
   emmap->derfa_z = NULL;
   emmap->dexpa_x = NULL;
   emmap->dexpa_y = NULL;
   emmap->dexpa_z = NULL;
   emmap->bound_x = NULL;
   emmap->bound_y = NULL;
   emmap->bound_z = NULL;
   emmap->erfa_x = NULL;
   emmap->erfa_y = NULL;
   emmap->erfa_z = NULL;
   emmap->expa_x = NULL;
   emmap->expa_y = NULL;
   emmap->expa_z = NULL;
   emmap->ig_x = NULL;
   emmap->ig_y = NULL;
   emmap->ig_z = NULL;
              
   read_emmap(emfile, emmap);
   vsize=emmap->nx*emmap->ny*emmap->nz;
//   snew(emmap->rho_em, vsize);
   snew(emmap->rho_sim, vsize);
   emmap->rdx=1/emmap->dx;
   emmap->rdy=1/emmap->dy;
   emmap->rdz=1/emmap->dz;   

   setCutoff(emmap);
   if(MASTER(cr)) {
     fprintf(stderr,"Cutoff size in x,y,z directions are %d %d %d\n",emmap->nGridCutX, emmap->nGridCutY,emmap->nGridCutZ);
     fprintf(fplog,"Cutoff size in x,y,z directions are %d %d %d\n",emmap->nGridCutX, emmap->nGridCutY,emmap->nGridCutZ);
   }
   
   snew(emmap->bound_x,emmap->nx+1);
   snew(emmap->bound_y,emmap->ny+1);
   snew(emmap->bound_z,emmap->nz+1);
   snew(emmap->erfa_x, emmap->nx+1);
   snew(emmap->erfa_y, emmap->ny+1);
   snew(emmap->erfa_z, emmap->nz+1);
   snew(emmap->expa_x, emmap->nx+1);
   snew(emmap->expa_y, emmap->ny+1);
   snew(emmap->expa_z, emmap->nz+1);
    
   expand_emmap(emmap, state,cr);

   for(i=0;i<=emmap->nx;i++) 
     emmap->bound_x[i]=emmap->x0 + emmap->dx*i - emmap->dx*0.5; 
   for(i=0;i<=emmap->ny;i++)
     emmap->bound_y[i]=emmap->y0 + emmap->dy*i - emmap->dy*0.5;
   for(i=0;i<=emmap->nz;i++)
     emmap->bound_z[i]=emmap->z0 + emmap->dz*i - emmap->dz*0.5;

   for(i=0;i<emmap->nx;i++) 
    for(j=0;j<emmap->ny;j++)
      for(k=0;k<emmap->nz;k++) {
         index = i* emmap->ny * emmap->nz + j *emmap->nz +k;
         norm_exp += emmap->rho_em[index] * emmap->rho_em[index];
         ave_exp += emmap->rho_em[index];
      }
   emmap->rho2_em = norm_exp;
   emmap->rhoave_em = ave_exp/(emmap->nx*emmap->ny*emmap->nz);
  if(MASTER(cr)){ 
     fprintf(stderr, "emmap initial rho^2 is %f\n", emmap->rho2_em);
     fprintf(fplog, "emmap initial rho^2 is %f\n", emmap->rho2_em);
  }
}

void free_emmap(t_emmap *emmap)
{
   sfree(emmap->derfa_x);
   sfree(emmap->derfa_y);
   sfree(emmap->derfa_z);
   sfree(emmap->dexpa_x);
   sfree(emmap->dexpa_y);
   sfree(emmap->dexpa_z);
   sfree(emmap->ig_x);
   sfree(emmap->ig_y);
   sfree(emmap->ig_z);
}

real update_sim_map(FILE *fplog,gmx_groups_t *group, t_state *state, rvec f[], t_emmap *emmap, t_mdatoms *md, t_commrec *cr, int mdfitgrp, int step)
{

//  printf("start update_sim_map\n");
   real *bound_x, *bound_y, *bound_z;
   real *erfa_x, *erfa_y, *erfa_z;
   real *derfa_x, *derfa_y, *derfa_z;
   real *expa_x, *expa_y, *expa_z;
   real *dexpa_x, *dexpa_y, *dexpa_z;
   int *ig_x, *ig_y, *ig_z;
   int n,i,j,k;
   int start, end, index, indey, indez, rhoidx,rhoidx2;
   real fcons1, fcons2;
   real drho_dx, drho_dy, drho_dz;
   real dot_exp_drhodx, dot_exp_drhody, dot_exp_drhodz;
   real dot_sim_drhodx, dot_sim_drhody, dot_sim_drhodz;
   real dcc_dx, dcc_dy, dcc_dz;
   int x1,x2,y1,y2,z1,z2;
   double norm_sim, ave_sim,dot_exp_sim,dot_exp_sim2,dot_exp_sim3,dot_exp_sim4;
   real cc,cc2,cc3,vcons;
   int flag, minx, maxx, miny, maxy, minz, maxz;
   
   start=0;
   end=state->natoms;
   
   /* the following allocation shall be done outside of this function*/
   bound_x=emmap->bound_x;
   bound_y=emmap->bound_y;
   bound_z=emmap->bound_z;
   erfa_x=emmap->erfa_x;
   erfa_y=emmap->erfa_y;
   erfa_z=emmap->erfa_z;
   derfa_x=emmap->derfa_x;
   derfa_y=emmap->derfa_y;
   derfa_z=emmap->derfa_z;
   expa_x=emmap->expa_x;
   expa_y=emmap->expa_y;
   expa_z=emmap->expa_z;
   dexpa_x=emmap->dexpa_x;
   dexpa_y=emmap->dexpa_y;
   dexpa_z=emmap->dexpa_z;
   ig_x=emmap->ig_x;
   ig_y=emmap->ig_y;
   ig_z=emmap->ig_z;
   emmap_zero(emmap); /* reset rho values to 0 memset(0)*/
   
   
   for(n=start; n<end; n++) {
     if (ggrpnr(group,egcMDFIT,n) < mdfitgrp) {
       ig_x[n]=floor( (state->x[n][0] - emmap->x0 + 0.5*emmap->dx)*emmap->rdx);
       ig_y[n]=floor( (state->x[n][1] - emmap->y0 + 0.5*emmap->dy)*emmap->rdy);
       ig_z[n]=floor( (state->x[n][2] - emmap->z0 + 0.5*emmap->dz)*emmap->rdz);
     }
   }
   
   ig_check(emmap, group, mdfitgrp, cr, state);
   
/* calculating partial rho */   
   flag = 0;

   for(n=start;n<end;n++) {
     if(ggrpnr(group,egcMDFIT,n) < mdfitgrp ) {
       x1 = ig_x[n] - emmap->nGridCutX;
       x2 = ig_x[n] + emmap->nGridCutX;
       y1 = ig_y[n] - emmap->nGridCutY;
       y2 = ig_y[n] + emmap->nGridCutY;
       z1 = ig_z[n] - emmap->nGridCutZ;
       z2 = ig_z[n] + emmap->nGridCutZ;

       if (flag == 0) {
         minx=x1;
         maxx=x2;
         miny=y1;
         maxy=y2;
         minz=z1;
         maxz=z2;
       }
       else {
         minx = (x1<minx ? x1:minx);
         maxx = (x2>maxx ? x2:maxx);
         miny = (y1<miny ? y1:miny);
         maxy = (y2>maxy ? y2:maxy);
         minz = (z1<minz ? z1:minz);
         maxz = (z2>maxz ? z2:maxz);
       }

       index=(n-start)*emmap->nx;
       indey=(n-start)*emmap->ny;
       indez=(n-start)*emmap->nz;
          
       for(j=x1;j<=x2;j++) {
         erfa_x[j] = erf(em_consErf*(bound_x[j] - state->x[n][0]));
         expa_x[j] = exp(em_consExp*(bound_x[j] - state->x[n][0])*(bound_x[j] - state->x[n][0]));
       }
       for(j=y1;j<=y2;j++) {
         erfa_y[j] = erf(em_consErf*(bound_y[j] - state->x[n][1]));
         expa_y[j] = exp(em_consExp*(bound_y[j] - state->x[n][1])*(bound_y[j] - state->x[n][1]));
       }       
       for(j=z1;j<=z2;j++) {
         erfa_z[j] = erf(em_consErf*(bound_z[j] - state->x[n][2]));
         expa_z[j] = exp(em_consExp*(bound_z[j] - state->x[n][2])*(bound_z[j] - state->x[n][2]));
       }       
       
       for(j=x1;j<x2;j++) {
         derfa_x[index+j]=erfa_x[j+1] - erfa_x[j];
         dexpa_x[index+j]=expa_x[j+1] - expa_x[j];
       }
       for(j=y1;j<y2;j++) {
         derfa_y[indey+j]=erfa_y[j+1] - erfa_y[j];
         dexpa_y[indey+j]=expa_y[j+1] - expa_y[j];
       }
       for(j=z1;j<z2;j++) {
         derfa_z[indez+j]=erfa_z[j+1] - erfa_z[j];
         dexpa_z[indez+j]=expa_z[j+1] - expa_z[j];
       }

       for(i=x1; i<x2; i++)
         for(j=y1; j<y2; j++)
           for(k=z1; k<z2; k++) {
	     rhoidx = i* emmap->ny * emmap->nz + j *emmap->nz +k;
	     emmap->rho_sim[rhoidx]+=derfa_x[index+i]*derfa_y[indey+j]*derfa_z[indez+k];
          }
    }
  }

   norm_sim=0.0;
   ave_sim=0.0;

   dot_exp_sim=0.0;
   dot_exp_sim2=0.0;
   dot_exp_sim3=0.0;
   dot_exp_sim4=0.0;

   
   vcons = pow(PI*emmap->sigma*emmap->sigma/6.0,3)/pow(emmap->dx*emmap->dy*emmap->dz,2);
   for(i=0;i<emmap->nx;i++)
     for(j=0;j<emmap->ny;j++)
       for(k=0;k<emmap->nz;k++) {
        rhoidx = i* emmap->ny * emmap->nz + j *emmap->nz +k;
	ave_sim+=emmap->rho_sim[rhoidx];
        norm_sim += emmap->rho_sim[rhoidx]*emmap->rho_sim[rhoidx];
	}

	ave_sim=ave_sim/(emmap->nx*emmap->ny*emmap->nz);

/* TRIPLE */    
   for(i=0;i<emmap->nx;i++)
     for(j=0;j<emmap->ny;j++)
       for(k=0;k<emmap->nz;k++) {
         rhoidx = i* emmap->ny * emmap->nz + j *emmap->nz +k;
	 dot_exp_sim += emmap->rho_em[rhoidx]*emmap->rho_sim[rhoidx];
	 dot_exp_sim2 += emmap->rho_em[rhoidx]*emmap->rho_sim[rhoidx];
	 dot_exp_sim3 += (emmap->rho_em[rhoidx]-emmap->rhoave_em)*(emmap->rho_sim[rhoidx] - ave_sim);
       }


   fcons1=1.0/sqrt(emmap->rho2_em * norm_sim * emmap->sigma);
   fcons2=dot_exp_sim/sqrt(emmap->rho2_em * pow(norm_sim,3) * emmap->sigma);


   for(n=start; n<end; n++) {
     if(ggrpnr(group,egcMDFIT,n) < mdfitgrp ) {
       dot_exp_drhodx = dot_exp_drhody = dot_exp_drhodz = 0;
       dot_sim_drhodx = dot_sim_drhody = dot_sim_drhodz = 0;

       x1 = ig_x[n] - emmap->nGridCutX;
       x2 = ig_x[n] + emmap->nGridCutX;
       y1 = ig_y[n] - emmap->nGridCutY;
       y2 = ig_y[n] + emmap->nGridCutY;
       z1 = ig_z[n] - emmap->nGridCutZ;
       z2 = ig_z[n] + emmap->nGridCutZ;
    /* TRIPLE */ 
       for (i = x1; i < x2; i++) {
         for (j = y1; j < y2; j++) {
           for (k = z1; k < z2; k++) {
	   
             index=(n-start)*emmap->nx;
             indey=(n-start)*emmap->ny;
             indez=(n-start)*emmap->nz;
	   
             drho_dx = dexpa_x[index+i] * derfa_y[indey+j] * derfa_z[indez+k];
	     drho_dy = dexpa_y[indey+j] * derfa_z[indez+k] * derfa_x[index+i];
	     drho_dz = dexpa_z[indez+k] * derfa_x[index+i] * derfa_y[indey+j];
	   
	     rhoidx = i* emmap->ny * emmap->nz + j *emmap->nz +k;	     
             dot_exp_drhodx += emmap->rho_em[rhoidx] * drho_dx; 
	     dot_sim_drhodx += emmap->rho_sim[rhoidx] * drho_dx;
	     dot_exp_drhody += emmap->rho_em[rhoidx] * drho_dy;
	     dot_sim_drhody += emmap->rho_sim[rhoidx] * drho_dy;
	     dot_exp_drhodz += emmap->rho_em[rhoidx] * drho_dz; 
	     dot_sim_drhodz += emmap->rho_sim[rhoidx] * drho_dz;
	  
	   }
         }
       }

       dcc_dx = dot_exp_drhodx * fcons1 - dot_sim_drhodx * fcons2;
       dcc_dy = dot_exp_drhody * fcons1 - dot_sim_drhody * fcons2;
       dcc_dz = dot_exp_drhodz * fcons1 - dot_sim_drhodz * fcons2;

       f[n][0] = emmap->weight * dcc_dx;
       f[n][1] = emmap->weight * dcc_dy;
       f[n][2] = emmap->weight * dcc_dz;
     }
   }

  cc  = dot_exp_sim/sqrt(emmap->rho2_em*norm_sim);
  cc3 = dot_exp_sim3/sqrt((emmap->rho2_em-emmap->rhoave_em*emmap->rhoave_em)*(norm_sim-ave_sim*ave_sim));
    if (MASTER(cr)){
      fprintf(stderr, "\nstep %d correlation coefficient: %f \n", step, cc);
      fprintf(fplog, "\nstep %d correlation coefficient: %f \n", step, cc);
}
  if(emmap->emwritefrequency >=1)
  if (step % emmap->emwritefrequency == 0)
    if (MASTER(cr)) {
      char filename[128];
      FILE *fp;
      sprintf(filename, "simmap.%d.sit",step);

      fprintf(stderr, "\nWriting simuated density at step %d to file %s \n", step, filename);
      fprintf(fplog, "\nWriting simuated density at step %d to file %s \n", step, filename);

      fp=fopen(filename,"w");
      fprintf(fp,"%10.6f %10.6f %10.6f %10.6f %8d %8d %8d\n\n",emmap->dx*10,emmap->x0*10,emmap->y0*10,emmap->z0*10,emmap->nx,emmap->ny,emmap->nz);
      int  stat = 0, vol;
/* TRIPLE */
      for(i=0;i<emmap->nx;i++)
        for(j=0;j<emmap->ny;j++) {
          for(k=0;k<emmap->nz;k++) {
             vol = k* emmap->ny * emmap->nz + j *emmap->nz +i;
             fprintf(fp,"%9.6f ",emmap->rho_sim[vol]);
	     stat++;
	     if (stat % 10 == 0)
	       fprintf(fp,"\n");
          }
        }    
      fclose(fp);
   } 
  return cc;

}

void update_sim_force(t_commrec *cr, rvec *em_f, rvec *f, t_mdatoms *md, int mdfitgrp)
{
  int start, end, i, ii;
  start=md->start;
  end=md->homenr + md->start;

  for(ii=0;ii<end;ii++) {
    i=cr->dd->gatindex[ii];
    if(md->cMDFIT[ii] < mdfitgrp) {
      f[ii][0] -=em_f[i][0];
      f[ii][1] -=em_f[i][1];
      f[ii][2] -=em_f[i][2];

    }
  }

}

void create_mdfit_comm(t_commrec *cr)
{
  int sim,i;
  int *rank;
#ifdef GMX_MPI
  MPI_Group mpi_group, mpi_n_group;

  cr->mpi_comm_mdfit = cr->mpi_comm_mysim;
#endif
  
  sim=cr->nnodes-1;
  snew(rank, sim);
  for(i=0;i<sim;i++) {
    rank[i]=i;
  }

#ifdef GMX_MPI  
  MPI_Comm_group(cr->mpi_comm_mdfit, &mpi_group);
  MPI_Group_incl(mpi_group, sim, rank, &mpi_n_group);
  MPI_Comm_create(cr->mpi_comm_mdfit, mpi_n_group, &cr->mpi_comm_mysim);
  cr->mpi_comm_mygroup = cr->mpi_comm_mysim;
#endif  
  sfree(rank);
}
