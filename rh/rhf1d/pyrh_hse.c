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

void _Hydrostatic(int NmaxIter, double iterLimit);

/* --- Global variables --                             -------------- */

Atmosphere atmos;
Geometry geometry;
Spectrum spectrum;
ProgramStats stats;
InputData input;
CommandLine commandline;
char messageStr[MAX_LINE_SIZE];

// enum Topology topology = ONE_D_PLANE;

/* --- Acceleration parameters --                      -------------- */

#define NG_HSE_DELAY   0
#define NG_HSE_ORDER   0
#define NG_HSE_PERIOD  0
#define NMAX_HSE_ITER  50

void dummy(){
  int argc = 1;
  char* argv[] = {"../rhf1d"};

  setOptions(argc, argv);
  SetFPEtraps();

  readInput();

  MULTIatmos(&atmos, &geometry);

  readAtomicModels();
}

/* ------- begin -------------------------- pyrh_hse.c -------------- */

myPops hse(int pyrh_Ndep, double pg_top,
           double *pyrh_scale, double *pyrh_temp, double *pyrh_ne, double *pyrh_vz, double *pyrh_vmic,
           double *pyrh_mag, double *pyrh_gamma, double *pyrh_chi,
           double *pyrh_nH, double *pyrh_nHtot, int pyrh_atm_scale, 
           int do_fudge, int fudge_num, double *fudge_lam, double *fudge)
{
  bool_t  equilibria_only, fromscratch;
  int     k, iter, index, layer;
  double  kappa, nHtot_old, eta, dcmass;
  double  muz, *S, *chi, *J, *rho, *pg, *chi_c;
  double  beta1, beta2;
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
  // Hydrogen populations are to be in LTE
  // It equals pointers of NLTE and LTE populations
  atmos.H_LTE = TRUE;

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
  // memcpy(atmos.ne, pyrh_ne, geometry.Ndep * sizeof(double));
  atmos.ne = malloc(atmos.Nspace * sizeof(double));
  memcpy(geometry.vel, pyrh_vz, geometry.Ndep * sizeof(double));
  memcpy(atmos.vturb, pyrh_vmic, geometry.Ndep * sizeof(double));
  memcpy(atmos.B, pyrh_mag, geometry.Ndep * sizeof(double));
  memcpy(atmos.gamma_B, pyrh_gamma, geometry.Ndep * sizeof(double));
  memcpy(atmos.chi_B, pyrh_chi, geometry.Ndep * sizeof(double));
  atmos.Stokes = TRUE;

  atmos.nH = matrix_double(atmos.NHydr, geometry.Ndep);
  
  // index=0;
  // for (int n=0; n<atmos.NHydr; n++)
  // {
  //   for (int k=0; k<geometry.Ndep; k++)
  //   {
  //     atmos.nH[n][k] = pyrh_nH[index];
  //     atmos.nH[n][k] /= CUBE(CM_TO_M);
  //     index++;
  //   }
  // }

  // set nHtot populations
  atmos.nHtot = (double *) calloc(geometry.Ndep, sizeof(double));
  
  // for (int k=0; k<geometry.Ndep; k++){
  //   atmos.nHtot[k]  = pyrh_nHtot[k];
  //   atmos.nHtot[k] /= CUBE(CM_TO_M);
  // }

  // check if atmosphere is non-static
  atmos.moving = FALSE;
  
  for (int k=0; k<geometry.Ndep; k++) {
    // for (int n=0;  n<atmos.NHydr;  n++) {
    //   atmos.nHtot[k] += atmos.nH[n][k];
    // }
    geometry.vel[k] *= KM_TO_M;
    atmos.vturb[k]  *= KM_TO_M;
    // atmos.ne[k]     /= CUBE(CM_TO_M);
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

  /* --- read atoms and molecules ----------------------------------- */
  
  readAtomicModels();
  readMolecularModels();
  
  /* --- set wavelength only to 500 nm ------------------------------ */
  
  double* wavetable = (double *) malloc(1 * sizeof(double));
  wavetable[0] = 500.00;
  int Nwav = 1;
  SortLambda(wavetable, Nwav);

  /* --- define variables for HSE------------------------------------ */

  myPops pops;
  pops.nH = matrix_double(6, atmos.Nspace);
  pops.ne = malloc(atmos.Nspace * sizeof(double));
  pops.nHtot = malloc(atmos.Nspace * sizeof(double));
  pops.rho = malloc(atmos.Nspace * sizeof(double));
  pops.pg = malloc(atmos.Nspace * sizeof(double));

  rho = (double *) malloc(atmos.Nspace * sizeof(double));
  pg = (double *) malloc(atmos.Nspace * sizeof(double));
  double* total_opacity = (double*) malloc(atmos.Nspace * sizeof(double)); // total opacity @ 500nm
  double* Nm            = (double*) malloc(atmos.Nspace * sizeof(double)); // total number density of molecules

  /*--- Start HSE solution for the top boundary  */

  // At the top boundary we specify gas pressure and obtain the electron pressure from it
  pg[0] = pg_top; // SI unit
  // printf("pyRH -- Pg = %e\n", pg[0]);
  
  atmos.active_layer = 0;
  iter = 0;

  // printf("%f | %e | %e \n", atmos.T[0], atmos.ne[0], atmos.nHtot[0]);

  nHtot_old = 1;
  atmos.nHmin[0] = 0;
  atmos.H2->n[0] = 0;
  atmos.ne[0] = 0;
  Nm[0] = 0;
  while (iter<NMAX_HSE_ITER){
    atmos.nHtot[0] = (pg[0]/KBOLTZMANN/atmos.T[0] - atmos.ne[0] - Nm[0]) / atmos.totalAbund + atmos.nHmin[0] + 2*atmos.H2->n[0];
    for (int n=0;  n<atmos.Natom; n++){
      atmos.atoms[n].ntotal[0] = atmos.atoms[n].abundance * atmos.nHtot[0];
    }
    rho[0] = (AMU * atmos.wght_per_H) * atmos.nHtot[0];
    get_ne(fromscratch=TRUE);
    SetLTEQuantities();
    pyrh_Background(equilibria_only=FALSE, total_opacity);
    get_Nm_total(Nm, atmos.active_layer);

    // pg[0] = atmos.gravity * geometry.tau_ref[0] * rho[0] / total_opacity[0];

    iter++;
    // eta = fabs((atmos.nHtot[0] - nHtot_old)/nHtot_old);
    // nHtot_old = atmos.nHtot[0];
    eta = fabs((atmos.ne[0] - nHtot_old)/nHtot_old);
    nHtot_old = atmos.ne[0];
    if (eta<=1e-2) break;
  }

  // if (eta>1e-2 && iter==NMAX_HSE_ITER) printf("pyRH -- Max number of iterations reached... eta = %e\n", eta);

  // printf("pyRH -- Pe = %e\n", atmos.ne[0] * KBOLTZMANN * atmos.T[0]*10);
  // printf("pyRH -- Pg = %e\n", pg[0]*10);

  // printf("%d | %f | %e | %e | %e | %e\n", iter, atmos.T[0], atmos.ne[0], atmos.nHtot[0], total_opacity[0], Nm[0]);

  double LOG10 = log(10);
  double dlogtau;

  /*--- Start HSE solution for rest atmospheric layers -------------- */

  for (k=1; k<atmos.Nspace; k++){
    iter = 0;
    atmos.active_layer = k;
    
    total_opacity[k] = total_opacity[k-1];
    rho[k] = rho[k-1];
    atmos.ne[k] = atmos.ne[k-1];
    Nm[k] = Nm[k-1];
    atmos.nHmin[k] = 0;
    atmos.H2->n[k] = 0;
    nHtot_old = 1;
    while (iter<NMAX_HSE_ITER){
      // de la Cruz Rodriguez et al. 2019
      // integration in logarithmic optical depth scale
      dlogtau = log10(geometry.tau_ref[k]) - log10(geometry.tau_ref[k-1]);
      if (iter==0){
        dcmass = 1 / (total_opacity[k]/rho[k]);
        // pg[k] = pg[k-1] + LOG10 * atmos.gravity * dlogtau * geometry.tau_ref[k] * dcmass;
        pg[k] = pg[k-1] + atmos.gravity * (geometry.tau_ref[k] - geometry.tau_ref[k-1]) * dcmass;
      }
      else{
        beta2 = total_opacity[k]/rho[k];
        beta1 = total_opacity[k-1]/rho[k-1];
        dcmass = log(beta2/beta1) / (beta2 - beta1);
        // dcmass = 2/(beta2 + beta1);
        // pg[k] = pg[k-1] + LOG10 * atmos.gravity * dlogtau * geometry.tau_ref[k] * dcmass;
        pg[k] = pg[k-1] + atmos.gravity * (geometry.tau_ref[k] - geometry.tau_ref[k-1]) * dcmass;
      }
      
      atmos.nHtot[k] = (pg[k]/KBOLTZMANN/atmos.T[k] - atmos.ne[k] - Nm[k]) / atmos.totalAbund + atmos.nHmin[k] + 2*atmos.H2->n[k];
      rho[k] = (AMU * atmos.wght_per_H) * atmos.nHtot[k];

      // get electron density and continuum opacity
      for (int n=0;  n<atmos.Natom; n++){
        atmos.atoms[n].ntotal[k] = atmos.atoms[n].abundance * atmos.nHtot[k];
      }
      get_ne(fromscratch=TRUE);
      SetLTEQuantities();
      pyrh_Background(equilibria_only=FALSE, total_opacity);
      get_Nm_total(Nm, atmos.active_layer);
      
      iter++;
      eta = fabs((atmos.nHtot[k] - nHtot_old)/nHtot_old);
      nHtot_old = atmos.nHtot[k];
      if (eta<=1e-2) break;

    }
    // if (eta>1e-2 && iter==NMAX_HSE_ITER) printf("pyRH -- Max number of iterations reached... eta = %e @ %d\n", eta, k);
    // printf("------------\n");
    // printf("%d | %e | %e | %e | %e | %e\n", iter, eta, atmos.ne[k], atmos.nHtot[k], total_opacity[k], Nm[k]);
    // printf("nHmin = %e\n", atmos.nHmin[k]);
    // printf("k = %d | iter = %d | eta = %f\n ----------- \n", k, iter, eta);
    // printf("k = %d | iter = %d\n ----------- \n", k, iter);
  }

  // store the populations
  pops.nH = atmos.atoms[0].nstar;
  pops.ne = atmos.ne;
  pops.nHtot = atmos.nHtot;
  pops.rho = memcpy(pops.rho, rho, atmos.Nspace * sizeof(double));
  pops.pg = memcpy(pops.pg, pg, atmos.Nspace * sizeof(double));
  
  // clear
  free(rho); rho = NULL;
  free(pg); pg = NULL;
  free(Nm); Nm = NULL;
  free(total_opacity); total_opacity = NULL;
  
  
  return pops;
}
/* ------- end ---------------------------- pyrh_hse.c -------------- */