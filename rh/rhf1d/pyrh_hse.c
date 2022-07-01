/* ------- file: -------------------------- solveray.c --------------

       Version:       rh2.0, 1-D plane-parallel
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Fri Jan 20 14:51:15 2012 --

       --------------------------                      ----------RH-- */

/* --- Solves radiative transfer for given atmosphere and model atom
       along a ray with arbitrary \mu_z, assuming the atom's population
       numbers and angle-averaged radiation field is given.


       Expects input file ``ray.input'' containing two lines of the form

         muz
         Nspect  wave_index1  ....   wave_indexNspect
       --                                              -------------- */

#include <fcntl.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "spectrum.h"
#include "statistics.h"
#include "inputs.h"
#include "error.h"
#include "xdr.h"

#include "constant.h"

#include "pyrh_hse.h"
#include "pyrh_background.h"

/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

// enum Topology topology = ONE_D_PLANE;

Atmosphere atmos;
Geometry geometry;
Spectrum spectrum;
ProgramStats stats;
InputData input;
CommandLine commandline;
char messageStr[MAX_LINE_SIZE];

void dummy(){
  int argc = 1;
  char* argv[] = {"../rhf1d"};

  setOptions(argc, argv);
  SetFPEtraps();

  readInput();

  MULTIatmos(&atmos, &geometry);

  readAtomicModels();
}

/* ------- begin -------------------------- solveray.c -------------- */

myPops hse(int pyrh_Ndep,
           double *pyrh_scale, double *pyrh_temp, double *pyrh_ne, double *pyrh_vz, double *pyrh_vmic,
           double *pyrh_mag, double *pyrh_gamma, double *pyrh_chi,
           double *pyrh_nH, int pyrh_atm_scale, 
           int do_fudge, int fudge_num, double *fudge_lam, double *fudge)
{
  bool_t  equilibria_only;
  int     k, iter, index, layer;
  double  kappa, nHtot_old, eta, dcmass;
  double  muz, *S, *chi, *J, *rho, *pg, *chi_c;
  Atom *atom;

  /* --- Read input data and initialize --             -------------- */
  int argc = 1;
  char* argv[] = {"../rhf1d"};

  setOptions(argc, argv);
  getCPU(0, TIME_START, NULL);
  SetFPEtraps();

  /* --- Read input data and initialize --             -------------- */

  readInput();
  spectrum.updateJ = TRUE;
  input.limit_memory = FALSE;
  // we want to solve for ne
  input.solve_ne = ONCE;
  input.startJ = NEW_J;
  
  // we want only to have reference wavelength
  // (but we will still have those from active atoms...)
  // input.wavetable_input = "none";
  // strcpy(input.wavetable_input, "none");
  // input.kurucz_list = "none";

  /* --- Read input data for atmosphere --             -------------- */

  if (input.StokesMode == FIELD_FREE ||
      input.StokesMode == POLARIZATION_FREE) {
    input.StokesMode = FULL_STOKES;
  }

  /* --- Setting up the atmosphere -- ------------------------------- */

  // set fudge factors
  if (do_fudge==1){
    input.do_fudge = TRUE;
    atmos.fudge_num = fudge_num;
    atmos.fudge_lam = fudge_lam;
    atmos.fudge = matrix_double(3, atmos.fudge_num);
    index = 0;
    for (int n=0; n<3; n++){
      for (int k=0; k<atmos.fudge_num; k++){
        atmos.fudge[n][k] = fudge[index];
        index++;
      }
    }
  }

  atmos.Nrlk = 0;

  geometry.Ndep = pyrh_Ndep;
  
  getCPU(1, TIME_START, NULL);
  MULTIatmos(&atmos, &geometry);
  
  if (pyrh_atm_scale==0){
    geometry.scale = TAU500;
    for (int k=0; k<geometry.Ndep; k++) 
      geometry.tau_ref[k] = POW10(pyrh_scale[k]);
  }
  if (pyrh_atm_scale==1){
    geometry.scale = COLUMN_MASS;
    for (int k=0; k<geometry.Ndep; k++){
      geometry.cmass[k] = POW10(pyrh_scale[k]) * (G_TO_KG / SQ(CM_TO_M));
    }
  }

  memcpy(atmos.T, pyrh_temp, geometry.Ndep * sizeof(double));
  memcpy(atmos.ne, pyrh_ne, geometry.Ndep * sizeof(double));
  memcpy(geometry.vel, pyrh_vz, geometry.Ndep * sizeof(double));
  memcpy(atmos.vturb, pyrh_vmic, geometry.Ndep * sizeof(double));
  memcpy(atmos.B, pyrh_mag, geometry.Ndep * sizeof(double));
  memcpy(atmos.gamma_B, pyrh_gamma, geometry.Ndep * sizeof(double));
  memcpy(atmos.chi_B, pyrh_chi, geometry.Ndep * sizeof(double));
  atmos.Stokes = TRUE;

  atmos.nH = matrix_double(atmos.NHydr, geometry.Ndep);
  index=0;
  for (int n=0; n<atmos.NHydr; n++)
  {
    for (int k=0; k<geometry.Ndep; k++)
    {
      atmos.nH[n][k] = pyrh_nH[index];
      atmos.nH[n][k] /= CUBE(CM_TO_M);
      index++;
    }
  }

  atmos.nHtot = (double *) calloc(geometry.Ndep, sizeof(double));
  
  // check if atmosphere is non-static
  atmos.moving = FALSE;
  
  for (int k=0; k<geometry.Ndep; k++) {
    for (int n=0;  n<atmos.NHydr;  n++) {
      atmos.nHtot[k] += atmos.nH[n][k];
    }
    geometry.vel[k] *= KM_TO_M;
    atmos.vturb[k]  *= KM_TO_M;
    atmos.ne[k]     /= CUBE(CM_TO_M);
  }

  for (int k=0; k<geometry.Ndep; k++)
  {
    if (fabs(geometry.vel[k]) >= atmos.vmacro_tresh) {
      atmos.moving = TRUE;
      break;
    }
  }

  /* --- redefine geometry for just this one ray --    -------------- */

  atmos.Nrays = geometry.Nrays = 1;
  geometry.muz[0] = 1.0;
  geometry.mux[0] = sqrt(1.0 - SQ(geometry.muz[0]));
  geometry.muy[0] = 0.0;
  geometry.wmu[0] = 1.0;
  if (atmos.Stokes) Bproject();

  readAtomicModels();
  readMolecularModels();
  double* wavetable = (double *) malloc(1 * sizeof(double));
  wavetable[0] = 500.00;
  int Nwav = 1;
  SortLambda(wavetable, Nwav);
  getBoundary(&geometry);

  myPops pops;
  pops.nH = matrix_double(6, atmos.Nspace);
  pops.ne = malloc(atmos.Nspace * sizeof(double));

  /*--- Start HSE solution for the top boundary  */

  rho = (double *) malloc(atmos.Nspace * sizeof(double));
  pg = (double *) malloc(atmos.Nspace * sizeof(double));
  double* total_opacity = (double*) malloc(atmos.Nspace * sizeof(double));

  iter = 0;
  atmos.active_layer = 0;
  // printf("nHtot = %e\n", atmos.nHtot[0]);
  while (iter<20){
    // get electron density and continuum opacity
    pyrh_Background(equilibria_only=FALSE, total_opacity);

    // get gas pressure
    rho[0] = (AMU * atmos.wght_per_H) * atmos.nHtot[0];
    kappa = total_opacity[0] / rho[0];
    dcmass = (geometry.tau_ref[0] / total_opacity[0]) * rho[0];
    pg[0] = atmos.gravity * dcmass;

    // convert it to nHtot
    nHtot_old = atmos.nHtot[0];
    atmos.nHtot[0] = (pg[0]/KBOLTZMANN/atmos.T[0] - atmos.ne[0]) / atmos.totalAbund;

    iter++;
    eta = fabs((atmos.nHtot[0] - nHtot_old)/nHtot_old);
    // printf(" eta = %e\n\n", eta);
    if (eta<=1e-2) break;
  }
  // printf("iter = %d | nHtot = %e\n", iter, atmos.nHtot[0]);

  // /*--- Start HSE solution for rest atmospheric layers  */

  for (k=1; k<atmos.Nspace; k++){
    // printf("ne = %e | %e | %e\n", atmos.ne[0], atmos.ne[1], atmos.ne[2]);
    iter = 0;
    atmos.active_layer = k;
    // printf("nHtot = %e\n", atmos.nHtot[k]);
    while (iter<20){
      // get electron density and continuum opacity
      pyrh_Background(equilibria_only=FALSE, total_opacity);

      // get gas pressure
      rho[k] = (AMU * atmos.wght_per_H) * atmos.nHtot[k];
      dcmass = (rho[k] + rho[k-1])/(total_opacity[k] + total_opacity[k-1]) * (geometry.tau_ref[k] - geometry.tau_ref[k-1]);
      pg[k] = pg[k-1] + atmos.gravity * dcmass;

      // convert it to nHtot
      nHtot_old = atmos.nHtot[k];
      atmos.nHtot[k] = (pg[k]/KBOLTZMANN/atmos.T[k] - atmos.ne[k]) / atmos.totalAbund;

      iter++;
      eta = fabs((atmos.nHtot[k] - nHtot_old)/nHtot_old);
      if (eta<=1e-2) break;
    }
    // printf("nHtot = %e\n", atmos.nHtot[k]);
    // printf("k = %d | iter = %d\n ----------- \n", k, iter);
  }
  // Final set-down of equilibrium values only
  pyrh_Background(equilibria_only=TRUE, total_opacity);

  /*--- Output the HSE atmosphere model ---*/

  free(rho); rho = NULL;
  free(pg); pg = NULL;
  free(total_opacity); total_opacity = NULL;
  
  pops.nH = atmos.atoms[0].nstar;
  pops.ne = atmos.ne;

  return pops;

  // printTotalCPU();

  // return;
}
/* ------- end ---------------------------- solveray.c -------------- */
