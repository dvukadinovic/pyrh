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

#include <stdlib.h>
// #include <stdio.h>
#include <string.h>

#include "../rh.h"
#include "../atom.h"
#include "../atmos.h"
#include "geometry.h"
#include "../spectrum.h"
#include "../background.h"
#include "../statistics.h"
#include "../error.h"
#include "../../headers/xdr.h"
#include "../constant.h"

#include "pyrh_compute1dray.h"

/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

enum Topology topology = ONE_D_PLANE;

/* --- Global variables --- */

Atmosphere atmos;
Geometry geometry;
Spectrum spectrum;
ProgramStats stats;
InputData input;
CommandLine commandline;
char messageStr[MAX_LINE_SIZE];

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

void concatenate(char* dest, char* str1, char* str2){
  strcpy(dest, str1);
  strcat(dest, str2);
}

myRLK_Line get_RLK_lines(char *cwd)
{
  myRLK_Line output;

  char* keyword_input = malloc(160);
  concatenate(keyword_input, cwd, "/keyword.input");
  strcpy(commandline.keyword_input, keyword_input);
  free(keyword_input);

  atmos.Stokes = TRUE;
  atmos.Nrlk = 0;

  setOptions(1, NULL);
  
  readInput();

  char* tmp = malloc(160);
  concatenate(tmp, "/", input.atoms_input);
  concatenate(input.atoms_input, cwd, tmp);
  free(tmp);

  // print out if we have ABO coeffs for each Kurucz line
  input.verbose = TRUE;

  readAbundance(&atmos, 0, NULL, NULL);
  // needs H atomic weight to compute ABO coeffs...
  readAtomicModels(); 

  readKuruczLines(input.KuruczData);

  output.rlk_lines = atmos.rlk_lines;
  output.Nrlk = atmos.Nrlk;

  return output;
}

mySpectrum rhf1d(char *cwd, double mu, int pyrh_Ndep,
              double *pyrh_scale, double *pyrh_temp, double *pyrh_ne, double *pyrh_vz, double *pyrh_vmic,
              double *pyrh_mag, double *pyrh_gamma, double *pyrh_chi,
              double *pyrh_nH, int pyrh_atm_scale, 
              int Nwave, double *lam,
              int fudge_num, double *fudge_lam, double *fudge,
              int Nloggf, int *loggf_ids, double *loggf_values,
              int Nlam, int *lam_ids, double *lam_values,
              int Nabun, int *atomic_id, double *atomic_abundance,
              int get_atomic_rfs, int get_populations,
              int NKurucz_lists, char *Kurucz_lists)
              // myRLK_Line *pyrh_rlk_lines,
{
  bool_t write_analyze_output, equilibria_only;
  int    niter, nact, index;
  Molecule *molecule;

  /* --- Read input data and initialize --             -------------- */
  char* keyword_input = malloc(160);
  concatenate(keyword_input, cwd, "/keyword.input");
  strcpy(commandline.keyword_input, keyword_input);
  free(keyword_input);

  setOptions(1, NULL);
  SetFPEtraps();

  readInput();

  input.verbose = FALSE;

  input.get_populations = get_populations;

  // We are not performing HSE; atoms and molecules can be NLTE
  input.pyrhHSE = FALSE;

  // get info if we need to compute semin-analytical RFs for atomic parameters
  input.get_atomic_rfs = FALSE;
  input.n_atomic_pars = 0;

  // input.solve_ne = ONCE;
  
  /*--- Overwrite values for ATOMS, MOLECULES and KURUCZ files ------ */
  char* tmp = malloc(160);
  
  // atomic list file
  concatenate(tmp, "/", input.atoms_input);
  concatenate(input.atoms_input, cwd, tmp);
  // molecules list file
  concatenate(tmp, "/", input.molecules_input);
  concatenate(input.molecules_input, cwd, tmp);
  // Kurucz list file
  concatenate(tmp, "/", input.KuruczData);
  concatenate(input.KuruczData, cwd, tmp);
  free(tmp);

  // DV: override/set some parameters 
  // atoms and molecules can be NLTE (if HSE=TRUE, than everything is LTE)
  input.pyrhHSE = FALSE;
  // we do not solve for ne pops, we provide them
  input.solve_ne = FALSE;
  // we have to recompute J
  spectrum.updateJ = TRUE;
  // do not print J in a file
  input.limit_memory = FALSE;
  // treatment of Hydrogen in the background. In LTE case it should be TRUE
  // atmos.H_LTE = TRUE; // --> move it after we check all atoms/molecules if they are active/passive?
  // if we need to solve NLTE problem
  input.solve_NLTE = FALSE;
  
  // --- DV --- this is where I stoped with reading the RH workflow.

  // set fudge factors
  if (fudge_lam!=NULL){
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
  if (Nloggf>=1){
    atmos.Nloggf = Nloggf;
    atmos.loggf_ids = loggf_ids;
    atmos.loggf_values = loggf_values;
    input.n_atomic_pars += Nloggf;
  }

  // set lam0 indices and values if forwarded
  if (Nlam>=1){
    atmos.Nlam = Nlam;
    atmos.lam_ids = lam_ids;
    atmos.lam_values = lam_values;
    input.n_atomic_pars += Nlam;
  }

  if (get_atomic_rfs!=0) input.get_atomic_rfs = TRUE;

  // if (pyrh_rlk_lines->Nrlk!=0){
  //   atmos.Nrlk = pyrh_rlk_lines->Nrlk;
  //   atmos.rlk_lines = pyrh_rlk_lines->rlk_lines;
  // } else {
  //   atmos.Nrlk = 0;
  // }
  atmos.Nrlk = 0;

  geometry.Ndep = pyrh_Ndep;
  
  getCPU(1, TIME_START, NULL);
  MULTIatmos(&atmos, &geometry, Nabun, atomic_id, atomic_abundance);
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
  if (pyrh_atm_scale==2){
    geometry.scale = GEOMETRIC;
    for (int k=0; k<geometry.Ndep; k++){
      geometry.height[k] = pyrh_scale[k] * KM_TO_M;
    }
  }

  atmos.T = pyrh_temp;
  atmos.ne = pyrh_ne;
  geometry.vel = pyrh_vz;
  atmos.vturb = pyrh_vmic;
  atmos.B = pyrh_mag;
  atmos.gamma_B = pyrh_gamma;
  atmos.chi_B = pyrh_chi;
  atmos.nHtot = pyrh_nH;
  atmos.nH = matrix_double(atmos.NHydr, geometry.Ndep);

  atmos.Stokes = TRUE;

  // check if atmosphere is non-static
  atmos.moving = FALSE;
  
  for (int k=0; k<geometry.Ndep; k++)
  {
    geometry.vel[k] *= KM_TO_M;
    atmos.vturb[k]  *= KM_TO_M;
    atmos.ne[k]     /= CUBE(CM_TO_M);
    atmos.nHtot[k]  /= CUBE(CM_TO_M);
    atmos.B[k] /= 1e4; // G --> T
  }

  for (int k=0; k<geometry.Ndep; k++)
  {
    if (fabs(geometry.vel[k]) >= atmos.vmacro_tresh) {
      atmos.moving = TRUE;
      break;
    }
  }

  readAtomicModels();
  readMolecularModels();

  // check if we can avoid NLTE computations
  // (2 calls for LTE solution are needed instead of only one)
  for (int n = 0;  n < atmos.Natom;  n++){
    if (atmos.atoms[n].active){
      input.solve_NLTE = TRUE;
      break;
    }
  }
  
  if (!input.solve_NLTE){
    atmos.Nrays = geometry.Nrays = 1;
    geometry.muz[0] = mu;
    geometry.mux[0] = sqrt(1.0 - SQ(geometry.muz[0]));
    geometry.muy[0] = 0.0;
    geometry.wmu[0] = 1.0;
  }
  if (atmos.Stokes) Bproject();

  SortLambda(lam, Nwave);

  // allocate space for atomic RFs if needed
  if (input.get_atomic_rfs){
    atmos.atomic_rfs = matrix3d_double(spectrum.Nspect, atmos.Nrays, input.n_atomic_pars);
  }

  getBoundary(&geometry);
  
  Background(write_analyze_output=FALSE, equilibria_only=FALSE);
  convertScales(&atmos, &geometry);
  // verifyed: pyrh and RH return the same tau scale from given populations (ne, nH)!

  // call in NLTE
  if (input.solve_NLTE) getProfiles();
  // here it initializes the spectrum and J;
  // it reads the J computed in previous run;
  // it also takes PRD data for lines from Active atoms;
  // watch it! J is connected with LIMIT_MEMORY keyword
  // call in NLTE
  bool_t pyrh_io_flag = FALSE;
  initSolution(pyrh_io_flag);
  if (input.solve_NLTE) initScatter();

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

  mySpectrum spec;
  spec.nlw = spectrum.Nspect;
  spec.Nrays = atmos.Nrays;

  _solveray(mu, &spec);

  free(spec.lam);
  free(spec.sI);
  free(spec.sQ);
  free(spec.sU);
  free(spec.sV);

  // revert units (since we pass pointers...)
  for (int k=0; k<geometry.Ndep; k++){
    geometry.vel[k] /= KM_TO_M;
    atmos.vturb[k]  /= KM_TO_M;
    atmos.ne[k]     *= CUBE(CM_TO_M);
    atmos.nHtot[k]  *= CUBE(CM_TO_M);
    atmos.B[k]      *= 1e4; // T --> G
  }

  //--- free all the memory that we do not use anymore

  if (atmos.moving || atmos.Stokes) free(atmos.backgrrecno);
  free(atmos.backgrflags);

  for (int idl=0; idl<atmos.Nrlk; idl++){
    free(atmos.rlk_lines[idl].zm->q);
    free(atmos.rlk_lines[idl].zm->strength);
    free(atmos.rlk_lines[idl].zm->shift);
    free(atmos.rlk_lines[idl].zm);
  }
  free(atmos.rlk_lines);
  
  freeAtoms();
  freeMolecules();
  freeElements();
  
  if (spectrum.wave_inds!=NULL) free(spectrum.wave_inds);
  if (spectrum.as!=NULL) free(spectrum.as); 

  if (atmos.Stokes){
    freeMatrix((void **) atmos.cos_gamma);
    freeMatrix((void **) atmos.cos_2chi);
    freeMatrix((void **) atmos.sin_2chi);
  }

  freeOpacityEmissivity();
  if (input.get_atomic_rfs) freeOpacityEmissivityDer();

  free(atmos.N);
  free(atmos.nHmin);

  // free geometry related data
  if (geometry.tau_ref!=NULL) free(geometry.tau_ref); geometry.tau_ref = NULL;
  if (geometry.cmass!=NULL) free(geometry.cmass); geometry.cmass = NULL;
  if (geometry.height!=NULL) free(geometry.height); geometry.height = NULL;

  if (geometry.Itop!=NULL) freeMatrix((void **) geometry.Itop);
  if (geometry.Ibottom!=NULL) freeMatrix((void **) geometry.Ibottom);

  free(geometry.wmu);
  free(geometry.mux);
  free(geometry.muy);
  free(geometry.muz);

  return spec;
}
/* ------- end ---------------------------- rhf1d.c ----------------- */
