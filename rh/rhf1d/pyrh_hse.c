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

#include "../rh.h"
#include "../atom.h"
#include "../atmos.h"
#include "geometry.h"
#include "../spectrum.h"
#include "../statistics.h"
#include "../inputs.h"
#include "../error.h"
#include "../../headers/xdr.h"

#include "../constant.h"

#include "pyrh_hse.h"
#include "pyrh_background.h"

/* --- Global variables --                             -------------- */

// Atmosphere atmos;
// Geometry geometry;
// Spectrum spectrum;
// InputData input;


/* --- Global variables --- */

extern Atmosphere atmos;
extern Geometry geometry;
extern Spectrum spectrum;
extern ProgramStats stats;
extern InputData input;
extern CommandLine commandline;
extern char messageStr[MAX_LINE_SIZE];

/* --- Acceleration parameters --                      -------------- */

#define NG_HSE_DELAY   0
#define NG_HSE_ORDER   0
#define NG_HSE_PERIOD  0
#define NMAX_HSE_ITER  50

/* ------- begin -------------------------- pyrh_hse.c -------------- */

void hse(char* cwd, int pyrh_Ndep,
           double *pyrh_scale, double *pyrh_temp,
           double *pyrh_ne, double *pyrh_nHtot, double *pyrh_rho, double *pyrh_pg,
           int pyrh_atm_scale, 
           int fudge_num, double *fudge_lam, double *fudge,
           int Nabun, int *abundance_id, double *abundance_value)
{
  bool_t  equilibria_only, fromscratch;
  int     k, iter, index, layer;
  double  kappa, nHtot_old, eta, dcmass;
  double  muz, *S, *chi, *J, *rho, *pg, *chi_c;
  double  beta1, beta2;
  Atom *atom;

  /* --- Read input data and initialize --             -------------- */
  // int argc = 1;
  // char* argv[] = {"../rhf1d"};//, "-i", keyword_input};
  // char* keyword_input = malloc(160);
  // concatenate(keyword_input, cwd, "/keyword.input");

  setOptions(1, NULL);
  SetFPEtraps();

  char* keyword_input = malloc(160);
  concatenate(keyword_input, cwd, "/keyword.input");
  strcpy(commandline.keyword_input, keyword_input);
  free(keyword_input);

  /* --- Read input data and initialize --             -------------- */

  readInput();
  // We are performing HSE; all atoms and molecules are to be treated in LTE
  input.pyrhHSE = TRUE;

  /*--- Overwrite values for ATOMS, MOLECULES and KURUCZ files ------ */
  char* tmp = malloc(160);
  
  // atomic list file
  concatenate(tmp, "/", input.atoms_input);
  concatenate(input.atoms_input, cwd, tmp);
  // molecules list file
  concatenate(tmp, "/", input.molecules_input);
  concatenate(input.molecules_input, cwd, tmp);
  free(tmp);

  spectrum.updateJ = TRUE;
  input.limit_memory = FALSE;
  // we want to solve for ne
  input.solve_ne = ONCE;
  input.startJ = NEW_J;
  // Hydrogen populations are to be in LTE
  // It equals pointers of NLTE and LTE populations
  atmos.H_LTE = TRUE;

  /* --- Read input data for atmosphere --             -------------- */

  if (input.StokesMode == FIELD_FREE ||
      input.StokesMode == POLARIZATION_FREE) {
    input.StokesMode = FULL_STOKES;
  }

  /* --- Setting up the atmosphere -- ------------------------------- */

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

  // no Kurucz lines
  atmos.Nrlk = 0;

  geometry.Ndep = pyrh_Ndep;

  getCPU(1, TIME_START, NULL);
  MULTIatmos(&atmos, &geometry, Nabun, abundance_id, abundance_value);
  
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
  atmos.nH = matrix_double(atmos.NHydr, geometry.Ndep);
  atmos.ne = pyrh_ne;
  atmos.nHtot = pyrh_nHtot;

  // it is needed in some routines...
  atmos.vturb = (double *) calloc(geometry.Ndep, sizeof(double));
  // atmos.vturb[0] = 2000; // [m/s]

  // no polarization
  atmos.Stokes = FALSE;
  // no velocities in the atmosphere
  atmos.moving = FALSE;

  /* --- redefine geometry for just this one ray --    -------------- */

  atmos.Nrays = geometry.Nrays = 1;
  geometry.muz[0] = 1.0;
  geometry.mux[0] = sqrt(1.0 - SQ(geometry.muz[0]));
  geometry.muy[0] = 0.0;
  geometry.wmu[0] = 1.0;

  /* --- read atoms and molecules ----------------------------------- */
  
  readAtomicModels();
  
  readMolecularModels();
  
  /* --- set wavelength only to 500 nm ------------------------------ */
  
  double* wavetable = (double *) malloc(1 * sizeof(double));
  wavetable[0] = 500.00;
  int Nwav = 1;
  SortLambda(wavetable, Nwav);
  
  /* --- define variables for HSE------------------------------------ */

  rho = pyrh_rho;
  pg = pyrh_pg;
  double* total_opacity = (double*) malloc(atmos.Nspace * sizeof(double)); // total opacity @ 500nm
  double* Nm            = (double*) malloc(atmos.Nspace * sizeof(double)); // total number density of molecules

  /*--- Start HSE solution for the top boundary  */

  atmos.backgrflags = NULL;
  atmos.active_layer = 0;
  iter = 0;

  // initial nH and ne
  atmos.ne[0] = 0;
  atmos.nHtot[0] = (pg[0]/KBOLTZMANN/atmos.T[0] - atmos.ne[0])/atmos.totalAbund;

  double deltaP;
  double turbP_el, turbP;

  // start iteration for the top
  while (iter<NMAX_HSE_ITER){
    // density
    rho[0] = (AMU * atmos.wght_per_H) * atmos.nHtot[0];

    // set ChemEq and LTE pops
    for (int n=0;  n<atmos.Natom; n++){
      atmos.atoms[n].ntotal[0] = atmos.atoms[n].abundance * atmos.nHtot[0];
    }
    get_ne(fromscratch=TRUE);
    SetLTEQuantities();
    
    // get opacity and molecule concentration
    pyrh_Background(equilibria_only=FALSE, total_opacity);
    get_Nm_total(Nm, atmos.active_layer);

    nHtot_old = atmos.nHtot[0];
    // atmos.nHtot[0] = (pg[0]/KBOLTZMANN/atmos.T[0] - atmos.ne[0] - Nm[0]) / atmos.totalAbund + atmos.nHmin[0] + 2*atmos.H2->n[0];
    // turbP_el = 0.5*M_ELECTRON*SQ(atmos.vturb[0])/KBOLTZMANN/atmos.T[0];
    // turbP = 0.5*atmos.avgMolWght*AMU*SQ(atmos.vturb[0])/KBOLTZMANN/atmos.T[0];
    // atmos.nHtot[0] = (pg[0]/KBOLTZMANN/atmos.T[0] - atmos.ne[0]*(1 + turbP_el)) / atmos.totalAbund / (1 + turbP);
    atmos.nHtot[0] = (pg[0]/KBOLTZMANN/atmos.T[0] - atmos.ne[0])/atmos.totalAbund;
  
    iter++;
    eta = fabs((atmos.nHtot[0] - nHtot_old)/atmos.nHtot[0]);
    if (eta<=1e-2) break;
  }
  
  if (geometry.scale==GEOMETRIC){
    geometry.cmass[0] = (atmos.nHtot[0] * atmos.totalAbund + atmos.ne[0]) * KBOLTZMANN*atmos.T[0] / atmos.gravity;
    geometry.tau_ref[0] = 0.5 * total_opacity[0] * (geometry.height[0] - geometry.height[1]);
    if (geometry.tau_ref[0] > 1.0) geometry.tau_ref[0] = 0.0;
  }

  double LOG10 = log(10);
  double dlogtau;

  /*--- Start HSE solution for the rest of atmospheric layers -------------- */

  for (k=1; k<atmos.Nspace; k++){
    iter = 0;
    atmos.active_layer = k;

    if (geometry.scale==GEOMETRIC){
      deltaP = atmos.gravity * rho[k-1] * (geometry.height[k-1] - geometry.height[k]);
    }
    if (geometry.scale==TAU500){
      deltaP = atmos.gravity * rho[k-1]/total_opacity[k-1] * (geometry.tau_ref[k] - geometry.tau_ref[k-1]);
    }
    pg[k] = pg[k-1] + deltaP;

    // turbP_el = 0.5*M_ELECTRON*SQ(atmos.vturb[k])/KBOLTZMANN/atmos.T[k];
    // turbP = 0.5*atmos.avgMolWght*AMU*SQ(atmos.vturb[k])/KBOLTZMANN/atmos.T[k];
    // atmos.nHtot[k] = (pg[k]/KBOLTZMANN/atmos.T[k] - atmos.ne[k]*(1 + turbP_el)) / atmos.totalAbund / (1 + turbP);
    atmos.nHtot[k] = (pg[k]/KBOLTZMANN/atmos.T[k] - atmos.ne[k-1])/atmos.totalAbund;

    while (iter<NMAX_HSE_ITER){
      rho[k] = (AMU * atmos.wght_per_H) * atmos.nHtot[k];
      
      // get electron density and LTE pops
      for (int n=0;  n<atmos.Natom; n++){
        atmos.atoms[n].ntotal[k] = atmos.atoms[n].abundance * atmos.nHtot[k];
      }
      get_ne(fromscratch=TRUE);
      SetLTEQuantities();

      // get continuum opacity and the molecular density
      pyrh_Background(equilibria_only=FALSE, total_opacity);
      get_Nm_total(Nm, atmos.active_layer);
      
      // de la Cruz Rodriguez et al. 2019
      // integration in logarithmic optical depth scale
      if (geometry.scale==TAU500) dlogtau = log10(geometry.tau_ref[k]) - log10(geometry.tau_ref[k-1]);
      
      if (geometry.scale==TAU500){
        beta1 = rho[k-1]/total_opacity[k-1] * geometry.tau_ref[k-1];
        beta2 = rho[k]/total_opacity[k] * geometry.tau_ref[k];
        pg[k] = pg[k-1] + LOG10 * atmos.gravity * (beta2 + beta1)/2 * dlogtau;
      }
      if (geometry.scale==GEOMETRIC){
        deltaP = atmos.gravity * (geometry.height[k-1] - geometry.height[k]) * sqrt(rho[k]*rho[k-1]);
        pg[k] = pg[k-1] + deltaP;
        
        // get other scales
        geometry.cmass[k]  = geometry.cmass[k-1]  + 0.5*(rho[k-1] + rho[k]) * (geometry.height[k-1] - geometry.height[k]);
        geometry.tau_ref[k] = geometry.tau_ref[k-1] + 0.5*(total_opacity[k-1] + total_opacity[k]) * (geometry.height[k-1] - geometry.height[k]);
      }

      nHtot_old = atmos.nHtot[k];
      // atmos.nHtot[k] = (pg[k]/KBOLTZMANN/atmos.T[k] - atmos.ne[k] - Nm[k]) / atmos.totalAbund + atmos.nHmin[k] + 2*atmos.H2->n[k];
      // turbP_el = 0.5*M_ELECTRON*SQ(atmos.vturb[k])/KBOLTZMANN/atmos.T[k];
      // turbP = 0.5*atmos.avgMolWght*AMU*SQ(atmos.vturb[k])/KBOLTZMANN/atmos.T[k];
      // atmos.nHtot[k] = (pg[k]/KBOLTZMANN/atmos.T[k] - atmos.ne[k]*(1 + turbP_el)) / atmos.totalAbund / (1 + turbP);
      atmos.nHtot[k] = (pg[k]/KBOLTZMANN/atmos.T[k] - atmos.ne[k]) / atmos.totalAbund;

      // printf("k = %d | %e | %e | %e\n", k, atmos.ne[k], atmos.nHtot[k], total_opacity[k]);
      iter++;
      eta = fabs((atmos.nHtot[k] - nHtot_old)/atmos.nHtot[k]);
      if (eta<=1e-2) break;
    }
  }

  // convertScales(&atmos, &geometry);

  //--- free all the memory that we do not use anymore

  free(atmos.backgrflags);

  if (geometry.Itop!=NULL) freeMatrix((void **) geometry.Itop);
  if (geometry.Ibottom!=NULL) freeMatrix((void **) geometry.Ibottom);

  free(Nm);
  free(total_opacity);

  ActiveSet *as;
  if (spectrum.wave_inds!=NULL) free(spectrum.wave_inds);
  // might have to free as at each wavelength first (for NLTE)
  for (int nspect = 0;  nspect < spectrum.Nspect;  nspect++) {
      as = &spectrum.as[nspect];
  
      free(as->Nlower);
      free(as->Nupper);

      free(as->upper_levels);
      free(as->lower_levels);
  }
  if (spectrum.as!=NULL) free(spectrum.as);

  free(wavetable);
  if (spectrum.lambda!=NULL) free(spectrum.lambda);
  
  freeOpacityEmissivity();
  if (input.get_atomic_rfs) freeOpacityEmissivityDer();

  if (atmos.Nactiveatom > 0) free(atmos.activeatoms);
  if (atmos.Nactivemol > 0) free(atmos.activemols);

  freeAtoms();
  freeMolecules();
  free(atmos.nHmin);
  
  free(atmos.vturb);
  free(atmos.N);
  
  free(geometry.wmu);
  free(geometry.mux);
  free(geometry.muy);
  free(geometry.muz);
  
  free(geometry.tau_ref);
  free(geometry.cmass);
  free(geometry.height);

  freeElements(); 
}
/* ------- end ---------------------------- pyrh_hse.c -------------- */

void get_scales(char *cwd, int pyrh_Ndep,
               double *pyrh_scale, double *pyrh_temp, double *pyrh_ne, double *pyrh_vz, double *pyrh_vmic,
               double *pyrh_nH, int pyrh_atm_scale, 
               double lam_ref, double *tau, double *height, double *cmass,
               int Nabun, int *abundance_id, double *abundance_value)
{
  /* --- Read input data and initialize --             -------------- */
  // int argc = 1;
  // char* argv[] = {"../rhf1d"};//, "-i", keyword_input};

  char* keyword_input = malloc(160);
  concatenate(keyword_input, cwd, "/keyword.input");
  strcpy(commandline.keyword_input, keyword_input);

  setOptions(1, NULL);
  SetFPEtraps();

  readInput();
  
  // We are performing 'HSE'; atoms and molecules are in LTE
  input.pyrhHSE = FALSE;
  
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

  input.StokesMode = NO_STOKES;

  atmos.Nrlk = 0;

  spectrum.updateJ = TRUE;
  input.startJ = NEW_J;
  input.limit_memory = FALSE;
  // For now, we only allow for H in LTE state
  atmos.H_LTE = TRUE;
  
  geometry.Ndep = pyrh_Ndep;
  
  getCPU(1, TIME_START, NULL);
  MULTIatmos(&atmos, &geometry, Nabun, abundance_id, abundance_value);

  if (pyrh_atm_scale==0){
    geometry.scale = TAU500;
    for (int k=0; k<geometry.Ndep; k++) geometry.tau_ref[k] = POW10(pyrh_scale[k]);
    geometry.height = height;
    geometry.cmass = cmass;
  }
  if (pyrh_atm_scale==1){
    geometry.scale = COLUMN_MASS;
    for (int k=0; k<geometry.Ndep; k++) geometry.cmass[k] = POW10(pyrh_scale[k]) * (G_TO_KG / SQ(CM_TO_M));
    geometry.height = height;
    geometry.tau_ref = tau;
  }
  if (pyrh_atm_scale==2){
    geometry.scale = GEOMETRIC;
    for (int k=0; k<geometry.Ndep; k++) geometry.height[k] = pyrh_scale[k] * KM_TO_M;
    geometry.cmass = cmass;
    geometry.tau_ref = tau;
  }

  atmos.T = pyrh_temp;
  atmos.ne = pyrh_ne;
  geometry.vel = pyrh_vz;
  atmos.vturb = pyrh_vmic;
  atmos.nHtot = pyrh_nH;
  atmos.nH = matrix_double(atmos.NHydr, geometry.Ndep);
 
  atmos.Stokes = FALSE;

  for (int k=0; k<geometry.Ndep; k++)
  {
    geometry.vel[k] *= KM_TO_M;
    atmos.vturb[k]  *= KM_TO_M;
    atmos.ne[k]     /= CUBE(CM_TO_M);
    atmos.nHtot[k]  /= CUBE(CM_TO_M);
  }

  // check if atmosphere is non-static
  atmos.moving = FALSE;
  for (int k=0; k<geometry.Ndep; k++)
  {
    if (fabs(geometry.vel[k]) >= atmos.vmacro_tresh) {
      atmos.moving = TRUE;
      break;
    }
  }
  
  readAtomicModels();
  readMolecularModels();

  double* wavetable = (double *) malloc(1 * sizeof(double));
  int Nwav = 1;
  wavetable[0] = lam_ref;
  atmos.lambda_ref = lam_ref;
  SortLambda(wavetable, Nwav);

  getBoundary(&geometry);

  atmos.active_layer = -1;
  input.solve_ne = NONE;

  Background(FALSE, FALSE);
  convertScales(&atmos, &geometry);

  // revert units (since we pass pointers...)
  for (int k=0; k<geometry.Ndep; k++){
    geometry.vel[k] /= KM_TO_M;
    atmos.vturb[k] /= KM_TO_M;
    atmos.ne[k] *= CUBE(CM_TO_M);
    atmos.nHtot[k] *= CUBE(CM_TO_M);
  }

  //--- free all the memory that we do not use anymore

  freeAtoms();
  freeMolecules();

  if (atmos.Stokes){
    freeMatrix((void **) atmos.cos_gamma);
    freeMatrix((void **) atmos.cos_2chi);
    freeMatrix((void **) atmos.sin_2chi);
  }

  freeOpacityEmissivity();
  freeElements();

  // free geometry related data
  // if (geometry.tau_ref!=NULL) free(geometry.tau_ref); geometry.tau_ref = NULL;
  // if (geometry.cmass!=NULL) free(geometry.cmass); geometry.cmass = NULL;
  // if (geometry.height!=NULL) free(geometry.height); geometry.height = NULL;

  if (spectrum.lambda!=NULL){
    free(spectrum.lambda); 
    spectrum.lambda = NULL;
    spectrum.Nspect = 0;
  }

  if (geometry.Itop!=NULL) freeMatrix((void **) geometry.Itop);
  if (geometry.Ibottom!=NULL) freeMatrix((void **) geometry.Ibottom);
}

void get_ne_from_nH(char *cwd, 
                    int pyrh_atm_scale, int pyrh_Ndep, 
                    double *pyrh_scale, double *pyrh_temp, 
                    double *pyrh_nH, double *pyrh_ne)
{
  /* --- Read input data and initialize --             -------------- */
  // int argc = 1;
  // char* argv[] = {"../rhf1d"};//, "-i", keyword_input};

  char* keyword_input = malloc(160);
  concatenate(keyword_input, cwd, "/keyword.input");
  strcpy(commandline.keyword_input, keyword_input);
  free(keyword_input);

  setOptions(1, NULL);
  SetFPEtraps();

  readInput();
  // input = pyrh_input;
  // We are not performing HSE; atoms and molecules can be NLTE
  input.pyrhHSE = TRUE;
  
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
  // input.KuruczData = NULL;

  free(tmp);

  atmos.Nrlk = 0;

  spectrum.updateJ = TRUE;
  input.limit_memory = FALSE;
  // For now, we only allow for H in LTE state
  atmos.H_LTE = TRUE;
  
  geometry.Ndep = pyrh_Ndep;
  
  getCPU(1, TIME_START, NULL);
  MULTIatmos(&atmos, &geometry, 0, NULL, NULL);

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
  atmos.nHtot = pyrh_nH;

  // 1/cm3 --> 1/m3
  for (int k=0; k<geometry.Ndep; k++){
    atmos.nHtot[k] /= CUBE(CM_TO_M);
  }

  // it is needed in some routines...
  atmos.vturb = (double *) calloc(geometry.Ndep, sizeof(double));

  atmos.Stokes = FALSE;
  atmos.moving = FALSE;

  atmos.nH = matrix_double(atmos.NHydr, geometry.Ndep);
  
  readAtomicModels();
  readMolecularModels();

  atmos.active_layer = -1;
  input.solve_ne = ONCE;
  Background(FALSE, TRUE);

  // 1/m3 --> 1/cm3
  for (int k=0; k<geometry.Ndep; k++){
    atmos.ne[k] *= CUBE(CM_TO_M);
    atmos.nHtot[k] *= CUBE(CM_TO_M);
  }

  //--- free all the memory that we do not use anymore

  freeAtoms();
  freeMolecules();

  if (atmos.Stokes){
    freeMatrix((void **) atmos.cos_gamma);
    freeMatrix((void **) atmos.cos_2chi);
    freeMatrix((void **) atmos.sin_2chi);
  }

  // free atmosphere related data
  free(atmos.vturb);
  freeElements();

  // free geometry related data
  if (geometry.tau_ref!=NULL) free(geometry.tau_ref); geometry.tau_ref = NULL;
  if (geometry.cmass!=NULL) free(geometry.cmass); geometry.cmass = NULL;
  if (geometry.height!=NULL) free(geometry.height); geometry.height = NULL;

  if (geometry.Itop!=NULL) freeMatrix((void **) geometry.Itop);
  if (geometry.Ibottom!=NULL) freeMatrix((void **) geometry.Ibottom);
}