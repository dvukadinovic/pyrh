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

#include "pyrh_compute1dray.h"
#include "pyrh_solveray.h"

#define COMMENT_CHAR        "#"

#define RAY_INPUT_FILE      "ray.input"
#define ASCII_SPECTRUM_FILE "spectrum_%4.2f.asc"

/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

// enum Topology topology = ONE_D_PLANE;

extern Atmosphere atmos;
extern Geometry geometry;
// ProgramStats stats;
extern InputData input;
// CommandLine commandline;

Spectrum spectrum;
char messageStr[MAX_LINE_SIZE];


int _getnumber(int* z)
{
  int a = 1;
  printf("z = %d\n", z[0]);
  return a+z[0];
}

/* ------- begin -------------------------- solveray.c -------------- */

// mySpectrum _solveray(int argc, char *argv[], int Ndep,
//               double *rh_scale, double *rh_temp, double *rh_ne, double *rh_vz, double *rh_vmic,
//               double *rh_mag, double *rh_gamma, double *rh_chi,
//               double **rh_nH, double muz,
//               int atm_scale)
void _solveray(char *argv[], double muz, mySpectrum *spec)
{
  register int n, k, la;

  char    rayFileName[14], inputLine[MAX_LINE_SIZE], ascFilename[18];
  bool_t  result, exit_on_EOF, to_obs, initialize, crosscoupling,
          analyze_output, equilibria_only;
  int     Nspect, Nread, Nrequired, checkPoint, *wave_index = NULL;
  double  *S, *chi, *J;
  FILE   *fp_out, *fp_ray, *fp_stokes, *fp_out_asc;
  XDR     xdrs;
  ActiveSet *as;

  // setOptions(argc, argv);
  // getCPU(0, TIME_START, NULL);
  // SetFPEtraps();

  /* --- Read input data and initialize --             -------------- */

  // readInput();
  input.startJ = OLD_J;
  spectrum.updateJ = FALSE;

  /* --- Read input data for atmosphere --             -------------- */

  // geometry.Ndep = Ndep;

  getCPU(1, TIME_START, NULL);
  // MULTIatmos(&atmos, &geometry);

  /* --- Read direction cosine for ray --              -------------- */

  // if ((fp_ray = fopen(RAY_INPUT_FILE, "r")) == NULL) {
  //   sprintf(messageStr, "Unable to open inputfile %s", RAY_INPUT_FILE);
  //   Error(ERROR_LEVEL_2, argv[0], messageStr);
  // }
  
  // getLine(fp_ray, COMMENT_CHAR, inputLine, exit_on_EOF=TRUE);
  // Nread = sscanf(inputLine, "%lf", &muz);
  // checkNread(Nread, Nrequired=1, argv[0], checkPoint=1);

  // if (muz <= 0.0  ||  muz > 1.0) {
  //   sprintf(messageStr,
	 //    "Value of muz = %f does not lie in interval <0.0, 1.0]\n", muz);
  //   Error(ERROR_LEVEL_2, argv[0], messageStr);
  // }
  
  // if (atm_scale==0){
  //   geometry.scale = TAU500;
  //   for (k=0; k<Ndep; k++){ 
  //     geometry.tau_ref[k] = POW10(rh_scale[k]);
  //   }
  // }
  // if (atm_scale==1){
  //   geometry.scale = COLUMN_MASS;
  //   for (k=0; k<Ndep; k++){
  //     geometry.cmass[k] = POW10(rh_scale[k]) * (G_TO_KG / SQ(CM_TO_M));
  //   }
  // }
  // if (input.StokesMode == FIELD_FREE ||
  //     input.StokesMode == POLARIZATION_FREE) {
  //   input.StokesMode = FULL_STOKES;
  // }
  /* --- redefine geometry for just this one ray --    -------------- */

  // atmos.T = rh_temp;
  // atmos.ne = rh_ne;
  // geometry.vel = rh_vz;
  // atmos.vturb = rh_vmic;

  // atmos.B = (double *) malloc(atmos.Nspace * sizeof(double));
  // atmos.gamma_B = (double *) malloc(atmos.Nspace * sizeof(double));
  // atmos.chi_B   = (double *) malloc(atmos.Nspace * sizeof(double));

  // atmos.B = rh_mag;
  // atmos.gamma_B = rh_gamma;
  // atmos.chi_B = rh_chi;
  // atmos.Stokes = TRUE;

  atmos.Nrays = geometry.Nrays = 1;
  geometry.muz[0] = muz;
  geometry.mux[0] = sqrt(1.0 - SQ(geometry.muz[0]));
  geometry.muy[0] = 0.0;
  geometry.wmu[0] = 1.0;
  // if (atmos.Stokes) Bproject();
  
  // atmos.nH = matrix_double(atmos.NHydr, Ndep);
  // atmos.nH = rh_nH;
  // atmos.nHtot = (double *) calloc(Ndep, sizeof(double));

  // for (k=0; k<Ndep; k++)
  // {
  //   for (n=0;  n<atmos.NHydr;  n++)
  //   {
  //     atmos.nH[n][k] /= CUBE(CM_TO_M);
  //     atmos.nHtot[k] += atmos.nH[n][k];
  //   }
  //   geometry.vel[k] *= KM_TO_M;
  //   atmos.vturb[k]  *= KM_TO_M;
  //   atmos.ne[k]     /= CUBE(CM_TO_M);
  // }
  
  // atmos.moving = FALSE;
  // for (k=0; k<Ndep; k++)
  // {
  //   if (fabs(geometry.vel[k]) >= atmos.vmacro_tresh) {
  //     atmos.moving = TRUE;
  //     break;
  //   }
  // }

  // this is how far I edited

  // readAtomicModels();
  // readMolecularModels();
  // SortLambda();

  // getBoundary(&geometry);

  /* --- Open file with background opacities --        -------------- */
  
  if (atmos.moving || input.StokesMode) {
    strcpy(input.background_File, input.background_ray_File);
    Background(analyze_output=FALSE, equilibria_only=FALSE);
  } else {
    Background(analyze_output=FALSE, equilibria_only=TRUE);

    if ((atmos.fd_background =
	 open(input.background_File, O_RDONLY, 0)) == -1) {
      sprintf(messageStr, "Unable to open inputfile %s",
	      input.background_File);
      Error(ERROR_LEVEL_2, argv[0], messageStr);
    }
    readBRS();
  }
  // convertScales(&atmos, &geometry);

  getProfiles();
  initSolution();
  initScatter();

  getCPU(1, TIME_POLL, "Total initialize");

  /* --- Solve radiative transfer equations --         -------------- */

  solveSpectrum(FALSE, FALSE);
  
  spec->nlw = spectrum.Nspect;
  spec->Nrays = atmos.Nrays;
  spec->lam = spectrum.lambda;
  spec->sI = spectrum.I[0];

  if (atmos.Stokes)
  {
    spec->sQ = spectrum.Stokes_Q[0];
    spec->sU = spectrum.Stokes_U[0];
    spec->sV = spectrum.Stokes_V[0];
    spec->stokes = 1;
  }
  else
  {
    spec->sQ = NULL;
    spec->sU = NULL;
    spec->sV = NULL;
    spec->stokes = 0;
  }
  spec->J = spectrum.J;
}

  /* --- Write emergent spectrum to output file --     -------------- */
 
  // sprintf(rayFileName, "spectrum_%4.2f", muz);
  // if ((fp_out = fopen(rayFileName, "w" )) == NULL) {
  //   sprintf(messageStr, "Unable to open output file %s", rayFileName);
  //   Error(ERROR_LEVEL_2, argv[0], messageStr);
  // }
  // xdrstdio_create(&xdrs, fp_out, XDR_ENCODE);

  // result = xdr_double(&xdrs, &muz);
  // result = xdr_vector(&xdrs, (char *) spectrum.I[0], spectrum.Nspect,
		//       sizeof(double), (xdrproc_t) xdr_double);

  /* --- Write ASCII table for special applications -- -------------- */

 //  if (!input.xdr_endian) {
 //    sprintf(ascFilename, ASCII_SPECTRUM_FILE, muz);
 //    if ((fp_out_asc = fopen(ascFilename, "w" )) == NULL) {
 //      sprintf(messageStr, "Unable to open output file %s", ascFilename);
 //      Error(ERROR_LEVEL_2, argv[0], messageStr);
 //    }
 //    fprintf(fp_out_asc, "%d\n", spectrum.Nspect);

 //    if (atmos.Stokes || input.backgr_pol) {
 //      for (la = 0;  la < spectrum.Nspect;  la++) {
	// fprintf(fp_out_asc, "%15.6lg %12.5lg %12.5lg %12.5lg %12.5lg\n",
	// 	spectrum.lambda[la],
	// 	spectrum.I[0][la], spectrum.Stokes_Q[0][la],
	// 	spectrum.Stokes_U[0][la], spectrum.Stokes_V[0][la]);
 //      }
 //    } else {
 //      for (la = 0;  la < spectrum.Nspect;  la++) {
	// fprintf(fp_out_asc, "%15.6lg %12.5lg %12.5lg %12.5lg %12.5lg\n",
	// 	spectrum.lambda[la],
	// 	spectrum.I[0][la], 0.0, 0.0, 0.0);
 //      }
 //    }
 //  }

  /* --- Read wavelength indices for which chi and S are to be
         written out for the specified direction --    -------------- */

  // Nread = fscanf(fp_ray, "%d", &Nspect);
  // checkNread(Nread, 1, argv[0], checkPoint=2);

  // if (Nspect > 0) {
  //   wave_index = (int *) malloc(Nspect * sizeof(int));
  //   Nread = 0;
  //   while (fscanf(fp_ray, "%d", &wave_index[Nread]) != EOF) Nread++;
  //   checkNread(Nread, Nspect, argv[0], checkPoint=3);
  //   fclose(fp_ray);

  //   chi = (double *) malloc(atmos.Nspace * sizeof(double));
  //   if (atmos.Stokes)
  //     S = (double *) malloc(4 * atmos.Nspace * sizeof(double));
  //   else
  //     S = (double *) malloc(atmos.Nspace * sizeof(double));
  // }
  // result = xdr_int(&xdrs, &Nspect);

  /* --- Go through the list of wavelengths --         -------------- */

  // if (Nspect > 0  &&  input.limit_memory)
  //   J = (double *) malloc(atmos.Nspace * sizeof(double));

  // for (n = 0;  n < Nspect;  n++) {
  //   if (wave_index[n] < 0  ||  wave_index[n] >= spectrum.Nspect) {
  //     sprintf(messageStr, "Illegal value of wave_index[n]: %4d\n"
	 //      "Value has to be between 0 and %4d\n", 
	 //      wave_index[n], spectrum.Nspect);
  //     Error(ERROR_LEVEL_2, argv[0], messageStr);
  //     continue;
  //   }
  //   sprintf(messageStr, "Processing n = %4d, lambda = %9.3f [nm]\n",
	 //    wave_index[n], spectrum.lambda[wave_index[n]]);
  //   Error(MESSAGE, NULL, messageStr);

  //   as = &spectrum.as[wave_index[n]];
  //   alloc_as(wave_index[n], crosscoupling=FALSE);
  //   Opacity(wave_index[n], 0, to_obs=TRUE, initialize=TRUE);
  //   readBackground(wave_index[n], 0, to_obs=TRUE);

  //   if (input.limit_memory) {
  //     readJlambda(wave_index[n], J);
  //   } else
  //     J = spectrum.J[wave_index[n]];

    /* --- Add the continuum opacity and emissivity -- -------------- */   

  //   for (k = 0;  k < atmos.Nspace;  k++) {
  //     chi[k] = as->chi[k] + as->chi_c[k];
  //     S[k]   = (as->eta[k] + as->eta_c[k] + as->sca_c[k]*J[k]) / chi[k];
  //   }
  //   result = xdr_int(&xdrs, &wave_index[n]);
  //   result = xdr_vector(&xdrs, (char *) chi, atmos.Nspace,
		// 	sizeof(double), (xdrproc_t) xdr_double);
  //   result = xdr_vector(&xdrs, (char *) S, atmos.Nspace,
		// 	sizeof(double), (xdrproc_t) xdr_double);

  //   free_as(wave_index[n], crosscoupling=FALSE);
  // }

  /* --- If magnetic fields are present --             -------------- */
  
  // if (atmos.Stokes || input.backgr_pol) {
  //   result = xdr_vector(&xdrs, (char *) spectrum.Stokes_Q[0],
		// 	spectrum.Nspect, sizeof(double),
		// 	(xdrproc_t) xdr_double);
  //   result = xdr_vector(&xdrs, (char *) spectrum.Stokes_U[0],
		// 	spectrum.Nspect, sizeof(double),
		// 	(xdrproc_t) xdr_double);
  //   result = xdr_vector(&xdrs, (char *) spectrum.Stokes_V[0],
		// 	spectrum.Nspect, sizeof(double),
		// 	(xdrproc_t) xdr_double);
  // }

  // if (Nspect > 0  &&  input.limit_memory)
  //   free(J);

  // xdr_destroy(&xdrs);
  // fclose(fp_out);

  // printTotalCPU();
// }
/* ------- end ---------------------------- solveray.c -------------- */
