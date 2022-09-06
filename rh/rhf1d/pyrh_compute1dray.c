/* ------- file: -------------------------- rhf1d.c -----------------

       Version:       rh2.0, 1-D plane-parallel
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Thu Feb 24 16:40:14 2011 --

       --------------------------                      ----------RH-- */

/* --- Main routine of 1D plane-parallel radiative transfer program.
       MALI scheme formulated according to Rybicki & Hummer

  See: G. B. Rybicki and D. G. Hummer 1991, A&A 245, p. 171-181
       G. B. Rybicki and D. G. Hummer 1992, A&A 263, p. 209-215

       Formal solution is performed with Feautrier difference scheme
       in static atmospheres, and with piecewise quadratic integration
       in moving atmospheres.

       --                                              -------------- */

// #include <stdlib.h>
// #include <stdio.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "spectrum.h"
#include "background.h"
#include "statistics.h"
#include "error.h"
#include "xdr.h"
#include "constant.h"

#include "pyrh_compute1dray.h"

/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

enum Topology topology = ONE_D_PLANE;

Atmosphere atmos;
Geometry geometry;
Spectrum spectrum;
ProgramStats stats;
InputData input;
CommandLine commandline;
char messageStr[MAX_MESSAGE_LENGTH];

/* ------- begin -------------------------- rhf1d.c ----------------- */

// We need in input here START_J key value:
//   if it is new, we then compute it
//   if it is old, we read it from the input
//
// The same is also true for the PRD, and atom/mol Populations
// (it they are started with OLD_POPULATIONS key; 
//  which is read from atoms.input file)
//
// These things should also be returned from rhf1d() in case we need it
// for later. Additionally, we should have separate calls to rhf1d() and 
// solveray(). Someone does not want to compute sovleray() maybe. If we 
// save every output from rhf1d() that we need for solveray(), then we can
// make it. For now, let us stick to solveray() call directly from the rhf1d().

myRLK_Line get_RLK_lines(int argc, char *argv[])
{
  myRLK_Line output;

  atmos.Stokes = TRUE;
  atmos.Nrlk = 0;
  setOptions(argc, argv);
  readInput();
  readAbundance(&atmos);
  readKuruczLines(input.KuruczData);

  output.rlk_lines = atmos.rlk_lines;
  output.Nrlk = atmos.Nrlk;

  return output;
}

// int argc, char *argv[], 
mySpectrum rhf1d(double mu, int pyrh_Ndep,
              double *pyrh_scale, double *pyrh_temp, double *pyrh_ne, double *pyrh_vz, double *pyrh_vmic,
              double *pyrh_mag, double *pyrh_gamma, double *pyrh_chi,
              double *pyrh_nH, int pyrh_atm_scale, 
              int Nwave, double *lam,
              int do_fudge, int fudge_num, double *fudge_lam, double *fudge,
              int Nloggf, int *loggf_ids, double *loggf_values,
              int Nlam, int *lam_ids, double *lam_values)
              // double *wavetable, int Nwave) // myRLK_Line *pyrh_rlk_lines,
{
  bool_t write_analyze_output, equilibria_only;
  int    niter, nact, index;

  Molecule *molecule;

  /* --- Read input data and initialize --             -------------- */
  int argc = 1;
  char* argv[] = {"../rhf1d"};

  setOptions(argc, argv);
  getCPU(0, TIME_START, NULL);
  SetFPEtraps();

  readInput();
  spectrum.updateJ = TRUE;
  input.limit_memory = FALSE;
  
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

  // set log(gf) indices and values if forwarded
  atmos.Nloggf = 0;
  atmos.loggf_ids = NULL;
  atmos.loggf_values = NULL;
  if (Nloggf>=1){
    atmos.Nloggf = Nloggf;
    atmos.loggf_ids = loggf_ids;
    atmos.loggf_values = loggf_values;
  }

  // set lam0 indices and values if forwarded
  atmos.Nlam = 0;
  atmos.lam_ids = NULL;
  atmos.lam_values = NULL;
  if (Nlam>=1){
    atmos.Nlam = Nlam;
    atmos.lam_ids = lam_ids;
    atmos.lam_values = lam_values;
  }

  // if (pyrh_rlk_lines->Nrlk!=0){
  //   atmos.Nrlk = pyrh_rlk_lines->Nrlk;
  //   atmos.rlk_lines = pyrh_rlk_lines->rlk_lines;
  // } else {
  //   atmos.Nrlk = 0;
  // }
  atmos.Nrlk = 0;

  geometry.Ndep = pyrh_Ndep;
  
  getCPU(1, TIME_START, NULL);
  MULTIatmos(&atmos, &geometry);
  atmos.active_layer = -1;
  
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
  
  for (int k=0; k<geometry.Ndep; k++)
  {
    for (int n=0;  n<atmos.NHydr;  n++)
    {
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

  if (atmos.Stokes) Bproject();
  
  readAtomicModels();
  readMolecularModels();
  SortLambda(lam, Nwave);

  getBoundary(&geometry);
  
  Background(write_analyze_output=FALSE, equilibria_only=FALSE);
  convertScales(&atmos, &geometry);

  bool_t pyrh_io_flag = FALSE;

  getProfiles();
  // here it initializes the spectrum 
  // and reads the J computed in rhf1d();
  // not only that, but takes PRD data also;
  // watch it! J is connected with LIMIT_MEMORY keyword
  initSolution(pyrh_io_flag);
  initScatter();

  getCPU(1, TIME_POLL, "Total Initialize");

  /* --- Solve radiative transfer for active ingredients -- --------- */

  // Here we get the spectrum (IQUV and J)
  Iterate(input.NmaxIter, input.iterLimit);

  adjustStokesMode();
  niter = 0;
  while (niter < input.NmaxScatter) {
    if (solveSpectrum(FALSE, FALSE) <= input.iterLimit) break;
    niter++;
  }

  if (atmos.hydrostatic) {
    geometry.scale = COLUMN_MASS;
    convertScales(&atmos, &geometry);
  }

  getCPU(1, TIME_START, NULL);

  mySpectrum spec;
  spec.nlw = Nwave;
  spec.Nrays = atmos.Nrays;

  // printf("J -- %e | %e \n", spectrum.J[0][0],  spectrum.J[0][10]);
  // printf("     %e | %e \n", spectrum.J[10][0], spectrum.J[10][10]);

  _solveray(argv, mu, &spec, spectrum.J, spectrum.J20);

  return spec;
}
/* ------- end ---------------------------- rhf1d.c ----------------- */
