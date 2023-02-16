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
#include <sys/sysinfo.h>

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

// Atmosphere atmos;
// Geometry geometry;
// Spectrum spectrum;
// ProgramStats stats;
// InputData input;
// CommandLine commandline;


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

InputData get_InputData(char *cwd)
{
  /* --- Read input data and initialize --             -------------- */
  int argc = 3;
  char* keyword_input = malloc(160);
  concatenate(keyword_input, cwd, "/keyword.input");
  char* argv[] = {"../rhf1d", "-i", keyword_input};

  setOptions(argc, argv);
  getCPU(0, TIME_START, NULL);
  SetFPEtraps();

  readInput();

  return input;
}

AtMol read_AtomsMolecules(InputData pyrh_input, char *cwd)
{
  AtMol AtomsMolecules;
  
  input = pyrh_input;
  
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

  readAbundance(&atmos);

  // AtomsMolecules.elements = atmos.elements;
  // AtomsMolecules.Nelem = atmos.Nelem;
  // AtomsMolecules.Tpf = atmos.Tpf;
  // AtomsMolecules.Npf = atmos.Npf;
  // AtomsMolecules.totalAbund = atmos.totalAbund;
  // AtomsMolecules.wght_per_H = atmos.avgWeight;
  // AtomsMolecules.avgMolWght = atmos.avgWeight / atmos.totalAbund;

  // readAtomicModels();
  // readMolecularModels();

  // AtomsMolecules.atoms = atmos.atoms;
  // AtomsMolecules.molecules = atmos.molecules;
  
  return AtomsMolecules;
}


void check_ID(InputData ID)
{
  printf("%d\n", ID.NmaxIter);
  printf("%s\n", ID.atoms_input);
}

// int argc, char *argv[], 
mySpectrum rhf1d(char *cwd, double mu, int pyrh_Ndep,
              double *pyrh_scale, double *pyrh_temp, double *pyrh_ne, double *pyrh_vz, double *pyrh_vmic,
              double *pyrh_mag, double *pyrh_gamma, double *pyrh_chi,
              double *pyrh_nH, int pyrh_atm_scale, 
              int Nwave, double *lam,
              int do_fudge, int fudge_num, double *fudge_lam, double *fudge,
              int Nloggf, int *loggf_ids, double *loggf_values,
              int Nlam, int *lam_ids, double *lam_values,
              int NKurucz_lists, char *Kurucz_lists)
              // myRLK_Line *pyrh_rlk_lines,
{
  bool_t write_analyze_output, equilibria_only;
  int    niter, nact, index;

  struct sysinfo sinfo;
  double totalRAM, freeRAM1, freeRAM2, mem_unit;
  int IOflag;

  // IOflag = sysinfo(&sinfo);
  // if (IOflag!=0){
  //   printf("Could not read sysinfo...\n");
  // } else{
  //   totalRAM = sinfo.totalram;
  //   mem_unit = sinfo.mem_unit / 1024.0 / 1024.0;
  //   freeRAM1 = sinfo.freeram * mem_unit;
  //   printf("free RAM %f\n", freeRAM1);
  // }

  Molecule *molecule;

  /* --- Read input data and initialize --             -------------- */
  int argc = 3;
  char* keyword_input = malloc(160);
  concatenate(keyword_input, cwd, "/keyword.input");
  char* argv[] = {"../rhf1d", "-i", keyword_input};

  // printf("NKurucz_lists = %d\n", NKurucz_lists);

  setOptions(argc, argv);
  getCPU(0, TIME_START, NULL);
  SetFPEtraps();

  readInput();
  // We are not performing HSE; atoms and molecules can be NLTE
  input.pyrhHSE = FALSE;
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

  spectrum.updateJ = TRUE;
  input.limit_memory = FALSE;
  // For now, we only allow for H in LTE state
  // atmos.H_LTE = TRUE;
  
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
  if (pyrh_atm_scale==2){
    geometry.scale = GEOMETRIC;
    for (int k=0; k<geometry.Ndep; k++){
      geometry.height[k] = pyrh_scale[k] * KM_TO_M;
    }
  }

  // memcpy(atmos.T, pyrh_temp, geometry.Ndep * sizeof(double));
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

  if (atmos.Stokes) Bproject();

  // IOflag = sysinfo(&sinfo);
  // if (IOflag!=0){
  //   printf("Could not read sysinfo...\n");
  // } else{
  //   freeRAM2 = sinfo.freeram * mem_unit;
  //   printf("After allocating atmosphere %f\n", freeRAM1 - freeRAM2);
  // }
  
  readAtomicModels();
  readMolecularModels();

  // IOflag = sysinfo(&sinfo);
  // if (IOflag!=0){
  //   printf("Could not read sysinfo...\n");
  // } else{
  //   freeRAM2 = sinfo.freeram * mem_unit;
  //   printf("Before allocating opacity/emisivity %f\n", freeRAM1 - freeRAM2);
  // }

  SortLambda(lam, Nwave);

  getBoundary(&geometry);

  // IOflag = sysinfo(&sinfo);
  // if (IOflag!=0){
  //   printf("Could not read sysinfo...\n");
  // } else{
  //   freeRAM2 = sinfo.freeram * mem_unit;
  //   printf("After allocating opacity/emisivity %f\n", freeRAM1 - freeRAM2);
  // }
  
  Background(write_analyze_output=FALSE, equilibria_only=FALSE);
  convertScales(&atmos, &geometry);
  // verifyed: pyrh and RH return the same tau scale from given populations (ne, nH)!

  // IOflag = sysinfo(&sinfo);
  // if (IOflag!=0){
  //   printf("Could not read sysinfo...\n");
  // } else{
  //   freeRAM2 = sinfo.freeram * mem_unit;
  //   printf("After Background() %f\n", freeRAM1 - freeRAM2);
  // }

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

  // IOflag = sysinfo(&sinfo);
  // if (IOflag!=0){
  //   printf("Could not read sysinfo...\n");
  // } else{
  //   freeRAM2 = sinfo.freeram * mem_unit;
  //   printf("Before iteration %f\n", freeRAM1 - freeRAM2);
  // }

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
  // spec.nlw = Nwave;
  spec.nlw = spectrum.Nspect;
  spec.Nrays = atmos.Nrays;

  // printf("J -- %e | %e \n", spectrum.J[0][0],  spectrum.J[0][10]);
  // printf("     %e | %e \n", spectrum.J[10][0], spectrum.J[10][10]);

  // IOflag = sysinfo(&sinfo);
  // if (IOflag!=0){
  //   printf("Could not read sysinfo...\n");
  // } else{
  //   freeRAM2 = sinfo.freeram * mem_unit;
  //   printf("Before solveray() %f\n", freeRAM1 - freeRAM2);
  // }

  _solveray(argv, mu, &spec);

  // revert units (since we pass pointers...)
  for (int k=0; k<geometry.Ndep; k++){
    geometry.vel[k] /= KM_TO_M;
    atmos.vturb[k]  /= KM_TO_M;
    atmos.ne[k]     *= CUBE(CM_TO_M);
    atmos.nHtot[k]  *= CUBE(CM_TO_M);
    atmos.B[k]      *= 1e4; // T --> G
  }

  //--- free all the memory that we do not use anymore
  for (int nspect=0; nspect < spectrum.Nspect; nspect++)
    free(spectrum.chi_c_lam[nspect]);
  free(spectrum.chi_c_lam); spectrum.chi_c_lam = NULL;
  
  if (input.magneto_optical){
    for (int nspect=0; nspect < spectrum.Nspect; nspect++)
      free(spectrum.chip_c_lam[nspect]);
    free(spectrum.chip_c_lam); spectrum.chip_c_lam = NULL;
  }

  for (int nspect=0; nspect < spectrum.Nspect; nspect++)
    free(spectrum.eta_c_lam[nspect]);
  free(spectrum.eta_c_lam); spectrum.eta_c_lam = NULL;

  for (int nspect=0; nspect < spectrum.Nspect; nspect++)
    free(spectrum.sca_c_lam[nspect]);
  free(spectrum.sca_c_lam); spectrum.sca_c_lam = NULL;

  // IOflag = sysinfo(&sinfo);
  // if (IOflag!=0){
  //   printf("Could not read sysinfo...\n");
  // } else{
  //   freeRAM2 = sinfo.freeram * mem_unit;
  //   printf("After freeing opacity/emisivity %f\n", freeRAM1 - freeRAM2);
  // }

  return spec;
}
/* ------- end ---------------------------- rhf1d.c ----------------- */
