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

InputData readMe(int argc, char *argv[])
{
  setOptions(argc, argv);
  readInput();

  return input;
}

mySpectrum rhf1d(int argc, char *argv[], int pyrh_Ndep,
              double *pyrh_scale, double *pyrh_temp, double *pyrh_ne, double *pyrh_vz, double *pyrh_vmic,
              double *pyrh_mag, double *pyrh_gamma, double *pyrh_chi,
              double **pyrh_nH, int pyrh_atm_scale)
{
  bool_t write_analyze_output, equilibria_only;
  int    niter, nact;

  Molecule *molecule;

  /* --- Read input data and initialize --             -------------- */
  setOptions(argc, argv);
  getCPU(0, TIME_START, NULL);
  SetFPEtraps();

  // readInput();
  spectrum.updateJ = TRUE;
  input.limit_memory = FALSE;

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
  
  // atmos.B = (double *) malloc(atmos.Nspace * sizeof(double));
  // atmos.gamma_B = (double *) malloc(atmos.Nspace * sizeof(double));
  // atmos.chi_B   = (double *) malloc(atmos.Nspace * sizeof(double));

  atmos.B = pyrh_mag;
  atmos.gamma_B = pyrh_gamma;
  atmos.chi_B = pyrh_chi;
  atmos.Stokes = TRUE;

  atmos.nH = matrix_double(atmos.NHydr, geometry.Ndep);
  for (int k=0; k<geometry.Ndep; k++)
  {
    for (int n=0; n<atmos.NHydr; n++)
    {
      atmos.nH[n][k] = pyrh_nH[n][k];
      atmos.nH[n][k] /= CUBE(CM_TO_M);
    }
  }
  atmos.nHtot = (double *) calloc(geometry.Ndep, sizeof(double));
  
  // check if atmosphere is non-static
  atmos.moving = FALSE;
  
  for (int k=0; k<geometry.Ndep; k++)
  {
    for (int n=0;  n<atmos.NHydr;  n++)
    {
      // atmos.nH[n][k] /= CUBE(CM_TO_M);
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
  SortLambda();
  
  getBoundary(&geometry);
  
  Background(write_analyze_output=TRUE, equilibria_only=FALSE);
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

  mySpectrum spec;
  _solveray(argv, 1.0, &spec, spectrum.J, spectrum.J20);

  return spec;
}
/* ------- end ---------------------------- rhf1d.c ----------------- */
