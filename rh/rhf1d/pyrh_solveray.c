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

// functions declaration
int _getnumber(int* z);
void _solveray(char *argv[], double muz, mySpectrum *spec, double** J, double** J20);


int _getnumber(int* z)
{
  int a = 1;
  printf("z = %d\n", z[0]);
  return a+z[0];
}

/* ------- begin -------------------------- solveray.c -------------- */

void _solveray(char *argv[], double muz, mySpectrum *spec, double** J, double** J20)
{
  register int n, k, la;

  char    rayFileName[14], inputLine[MAX_LINE_SIZE], ascFilename[18];
  bool_t  result, exit_on_EOF, to_obs, initialize, crosscoupling,
          analyze_output, equilibria_only;
  int     Nspect, Nread, Nrequired, checkPoint, *wave_index = NULL;
  double  *S, *chi;
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
  input.limit_memory = FALSE;

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


  bool_t pyrh_io_flag = FALSE;

  getProfiles();
  initSolution(pyrh_io_flag);
  spectrum.J = J;
  if (input.backgr_pol) spectrum.J20 = J20;
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