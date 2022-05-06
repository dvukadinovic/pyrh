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
mySpectrum rhf1d(int pyrh_Ndep,
              double *pyrh_scale, double *pyrh_temp, double *pyrh_ne, double *pyrh_vz, double *pyrh_vmic,
              double *pyrh_mag, double *pyrh_gamma, double *pyrh_chi,
              double *pyrh_nH, int pyrh_atm_scale,
              double *wavetable, int Nwave) // myRLK_Line *pyrh_rlk_lines,
{
  bool_t write_analyze_output, equilibria_only;
  int    niter, nact;

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

  // if (pyrh_rlk_lines->Nrlk!=0){
  //   atmos.Nrlk = pyrh_rlk_lines->Nrlk;
  //   atmos.rlk_lines = pyrh_rlk_lines->rlk_lines;
  // } else {
  //   atmos.Nrlk = 0;
  // }
  atmos.Nrlk = 0;

  Atom *atom;

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

  atmos.T = pyrh_temp;
  atmos.ne = pyrh_ne;
  geometry.vel = pyrh_vz;
  atmos.vturb = pyrh_vmic;
  atmos.B = pyrh_mag;
  atmos.gamma_B = pyrh_gamma;
  atmos.chi_B = pyrh_chi;
  atmos.Stokes = TRUE;

  atmos.nH = matrix_double(atmos.NHydr, geometry.Ndep);
  int index=0;
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
  SortLambda(wavetable, Nwave);

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
  /* --- Write output files --                         -------------- */

  if (atmos.hydrostatic) {
    geometry.scale = COLUMN_MASS;
    convertScales(&atmos, &geometry);
  }

  getCPU(1, TIME_START, NULL);

  // writeInput();
  // writeAtmos(&atmos);
  // writeGeometry(&geometry);
  // writeSpectrum(&spectrum);
  // writeFlux(FLUX_DOT_OUT);

  // for (nact = 0;  nact < atmos.Nactiveatom;  nact++) {
  //   atom = atmos.activeatoms[nact];

  //   writeAtom(atom);
  //   writePopulations(atom);
  //   writeRadRate(atom);
  //   writeCollisionRate(atom);
  //   writeDamping(atom);
  // }

  // for (nact = 0;  nact < atmos.Nactivemol;  nact++) {
  //   molecule = atmos.activemols[nact];
  //   writeMolPops(molecule);
  // }

  // writeOpacity();
  
  // getCPU(1, TIME_POLL, "Write output");

  printTotalCPU();

  /*--- Free from memory background opacity/emissivity  ---*/
  
  // int nspect;
  
  // for (nspect=0; nspect<spectrum.Nspect; nspect++) {
  //   free(spectrum.chi_c_lam[nspect]);
  // }
  // free(spectrum.chi_c_lam);

  // for (int nspect=0; nspect<spectrum.Nspect; ++nspect) {
  //   free(spectrum.eta_c_lam[nspect]);
  // }
  // free(spectrum.eta_c_lam);

  // for (int nspect=0; nspect<spectrum.Nspect; ++nspect) {
  //   free(spectrum.sca_c_lam[nspect]);
  // }
  // free(spectrum.sca_c_lam);

  // if (input.magneto_optical){
  //   for (int nspect=0; nspect<spectrum.Nspect; ++nspect) {
  //     free(spectrum.chip_c_lam[nspect]);
  //   }
  //   free(spectrum.chip_c_lam);
  // }

  mySpectrum spec;
  _solveray(argv, 1.0, &spec, spectrum.J, spectrum.J20);

  return spec;
}
/* ------- end ---------------------------- rhf1d.c ----------------- */
