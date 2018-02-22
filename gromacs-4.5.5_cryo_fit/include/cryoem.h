#ifndef _CRYOEM_H
#define _CRYOEM_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

//#ifdef GMX_MPI
//#include <mpi.h>
//#endif

#include "typedefs.h"

// Stuff for MDfit

typedef struct {
  int nx, ny, nz; /* number of grids in x,y,z direction*/
  int nxV, nyV, nzV; /* number of grids in x,y,z direction*/
  real dx, dy, dz; /* size of grid */
  real dxV, dyV, dzV; /* size of grid */
  real rdx, rdy, rdz; /* inverse of dx, dy, dz */
  real x0, y0, z0; /* center of the grids */
  real x0V, y0V, z0V; /* center of the grids */
  real *rho_em, *rho_sim ; /*density rho[nx][ny][nz] */
  real sigma, err, weight; /* Gaussian sigma -resolution; cutoff distance; K value(Potential weight) */
  int steps; /* frequency of using MDFIT correction */
  double rho2_em; /* sum_ijk(rho*rho) */
  double rhoave_em; /* sum_ijk(rho*rho) */
  int nGridCutX, nGridCutY, nGridCutZ; /* number of grids within cutoff distance */
  /* the following are the auxiliary data structures to calculate rho and rho' */
  real *bound_x, *bound_y, *bound_z;  /* size: nx*ny*nz */
  real *erfa_x, *erfa_y, *erfa_z; /* size: nx+1, ny+1, nz+1 */
  real *derfa_x, *derfa_y, *derfa_z; /* size: natoms*nx, natoms*ny, natoms*nz */
  real *expa_x, *expa_y, *expa_z; /* size: nx+1, ny+1, nz+1 */
  real *dexpa_x, *dexpa_y, *dexpa_z; /* size: nr*nx, nr*ny, nr*nz */
  int *ig_x, *ig_y, *ig_z; /* size: natoms */
  int nr; /*number of local atoms */
  real thresh; /*number of local atoms */
  int emwritefrequency; /* frequency for writing out theoretical maps*/ 
} t_emmap;

/* initialization of emmap data structure */
extern void init_em_map(FILE *fplog, char *emfile,t_emmap *emmap, t_inputrec *ir, t_state *state, t_commrec *cr);

/* reallocate space */
extern void expand_emmap(t_emmap *emmap, t_state *state,  t_commrec *cr);

/* free any allocated memory for emmap */
extern void free_emmap(t_emmap *emmap);

/* update force by MDFIT potential */
extern real update_sim_map(FILE *fplog,gmx_groups_t *group, t_state *state, rvec f[], t_emmap *emmap, t_mdatoms *md, t_commrec *cr,int mdfitgrp, int step);

extern void update_sim_force(t_commrec *cr, rvec *em_f, rvec *f, t_mdatoms *md, int mdfitgrp);

#endif
