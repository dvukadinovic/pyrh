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
#include "background.h"
#include "statistics.h"
#include "inputs.h"
#include "error.h"
#include "xdr.h"

#include "constant.h"

#include "pyrh_hse.h"

#define COMMENT_CHAR        "#"

#define RAY_INPUT_FILE      "ray.input"
#define ASCII_SPECTRUM_FILE "spectrum_%4.2f.asc"


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

/* ------- begin -------------------------- solveray.c -------------- */

myPops hse(int argc, char *argv[])
{
  register int n, k, la;

  char    rayFileName[14], inputLine[MAX_LINE_SIZE], ascFilename[18];
  bool_t  result, exit_on_EOF, to_obs, initialize, crosscoupling;
  bool_t  write_analyze_output, equilibria_only;
  int     Nspect, Nread, Nrequired, checkPoint, *wave_index = NULL, ref_index, iter;
  double  kappa, nHtot_old, eta, dcmass;
  double  muz, *S, *chi, *J, *rho, *pg;
  FILE   *fp_out, *fp_ray, *fp_stokes, *fp_out_asc;
  ActiveSet *as;

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

  /* --- Read input data for atmosphere --             -------------- */

  getCPU(1, TIME_START, NULL);
  oldMULTIatmos(&atmos, &geometry);
  rho = (double *) malloc(atmos.Nspace * sizeof(double));
  pg = (double *) malloc(atmos.Nspace * sizeof(double));

 if (input.StokesMode == FIELD_FREE ||
      input.StokesMode == POLARIZATION_FREE) {
    input.StokesMode = FULL_STOKES;
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
  SortLambda(NULL, 0);
  getBoundary(&geometry);

  // convertScales(&atmos, &geometry);

  // pops->nH = matrix_double(atmos.Nspace, atmos.Nspace);
  // pops->ne = malloc(atmos.Nspace * sizeof(double));

  /*--- Start HSE solution for the top boundary  */

  iter = 0;
  while (iter<20){
    // get electron density and continuum opacity
    Background(write_analyze_output=FALSE, equilibria_only=FALSE, 0);

    // get opacity at reference wavelength (500nm)
    Locate(spectrum.Nspect, spectrum.lambda, atmos.lambda_ref, &ref_index);

    as = &spectrum.as[ref_index];
    alloc_as(ref_index, FALSE);
    readBackground(ref_index, 0, 0);

    // get gas pressure
    rho[0] = (AMU * atmos.wght_per_H) * atmos.nHtot[0];
    kappa = as->chi_c[0] / rho[0];
    dcmass = (geometry.tau_ref[0] / as->chi_c[0]) * rho[0];
    pg[0] = atmos.gravity * dcmass;

    // convert it to nHtot
    nHtot_old = atmos.nHtot[0];
    atmos.nHtot[0] = (pg[0]/KBOLTZMANN/atmos.T[0] - atmos.ne[0]) / atmos.totalAbund;

    iter++;
    eta = fabs((atmos.nHtot[0] - nHtot_old)/nHtot_old);
    if (eta<=1e-2) break;
  }

  /*--- Start HSE solution for rest atmospheric layers  */

  for (k=1; k<atmos.Nspace; k++){
    iter = 0;
    while (iter<20){
      // get electron density and continuum opacity
      Background(write_analyze_output=FALSE, equilibria_only=FALSE, k);

      // get opacity at reference wavelength (500nm)
      // Locate(spectrum.Nspect, spectrum.lambda, atmos.lambda_ref, &ref_index);

      as = &spectrum.as[ref_index];
      alloc_as(ref_index, FALSE);
      readBackground(ref_index, 0, 0);

      // get gas pressure
      rho[k] = (AMU * atmos.wght_per_H) * atmos.nHtot[k];
      dcmass = (rho[k] + rho[k-1])/(as->chi_c[k] + as->chi_c[k-1]) * (geometry.tau_ref[k] - geometry.tau_ref[k-1]);
      pg[k] = pg[k-1] + atmos.gravity * dcmass;

      // convert it to nHtot
      nHtot_old = atmos.nHtot[k];
      atmos.nHtot[k] = (pg[k]/KBOLTZMANN/atmos.T[k] - atmos.ne[k]) / atmos.totalAbund;

      iter++;
      eta = fabs((atmos.nHtot[k] - nHtot_old)/nHtot_old);
      if (eta<=1e-2) break;
    }
    // printf("k = %d | iter = %d\n", k, iter);
    // printf("nHtot = %e\n", atmos.nHtot[k]);
  }

  /*--- Output the HSE atmosphere model ---*/

  // free(rho);
  // free(pg);

  myPops pops;
  pops.nH = atmos.H->n;
  pops.ne = atmos.ne;

  return pops;

  // printTotalCPU();
}
/* ------- end ---------------------------- solveray.c -------------- */
