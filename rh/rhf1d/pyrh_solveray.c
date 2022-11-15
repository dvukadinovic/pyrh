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

  /* --- Read input data and initialize --             -------------- */

  input.startJ = OLD_J;
  spectrum.updateJ = FALSE;
  input.limit_memory = FALSE;

  /* --- Read input data for atmosphere --             -------------- */

  getCPU(1, TIME_START, NULL);
  
  atmos.Nrays = geometry.Nrays = 1;
  geometry.muz[0] = muz;
  geometry.mux[0] = sqrt(1.0 - SQ(geometry.muz[0]));
  geometry.muy[0] = 0.0;
  geometry.wmu[0] = 1.0;
  // this has to be reprojected for new muz (but from input B vector, not from already projected from rhf1d())
  if (atmos.Stokes) Bproject();

  /* --- Open file with background opacities --        -------------- */

  if (atmos.moving || input.StokesMode) {
    // strcpy(input.background_File, input.background_ray_File);
    Background(analyze_output=FALSE, equilibria_only=FALSE);
  } else {
    Background(analyze_output=FALSE, equilibria_only=TRUE);

  //   if ((atmos.fd_background =
	 // open(input.background_File, O_RDONLY, 0)) == -1) {
  //     sprintf(messageStr, "Unable to open inputfile %s",
	 //      input.background_File);
  //     Error(ERROR_LEVEL_2, argv[0], messageStr);
  //   }
  //   readBRS();
  }
  // convertScales(&atmos, &geometry);


  bool_t pyrh_io_flag = FALSE;

  getProfiles();
  // spectrum.J is already filled with correct values; do we need this initSolution() here?
  // initSolution(pyrh_io_flag);
  // spectrum.J = J;
  // if (input.backgr_pol) spectrum.J20 = J20;
  initScatter();

  getCPU(1, TIME_POLL, "Total initialize");

  /* --- Solve radiative transfer equations --         -------------- */

  solveSpectrum(FALSE, FALSE);
  
  spec->sQ = NULL;
  spec->sU = NULL;
  spec->sV = NULL;
  spec->J = NULL;
  spec->J20 = NULL;

  int Nlw = spec->nlw;// - 1;
  Nlw -= 1;

  spec->stokes = 0;
  spec->lam = (double *) malloc(Nlw * sizeof(double));
  spec->sI = (double *) malloc(Nlw * sizeof(double));
  spec->sQ = (double *) malloc(Nlw * sizeof(double));
  spec->sU = (double *) malloc(Nlw * sizeof(double));
  spec->sV = (double *) malloc(Nlw * sizeof(double));
  // spec->J  = matrix_double(Nlw+1, atmos.Nspace);

  int index=0;
  double tmp;
  
  for (int idl=0; idl<Nlw+1; idl++){
    // skip referent wavelength
    if (spectrum.lambda[idl]==atmos.lambda_ref) continue;
    vacuum_to_air(1, &spectrum.lambda[idl], &tmp);
    spec->lam[index] = tmp;
    spec->sI[index] = spectrum.I[idl][0];
    if (atmos.Stokes){
      spec->sQ[index] = spectrum.Stokes_Q[idl][0];
      spec->sU[index] = spectrum.Stokes_U[idl][0];
      spec->sV[index] = spectrum.Stokes_V[idl][0];
      spec->stokes = 1;
    }
    index += 1;
    // for (int idz=0; idz<atmos.Nspace; idz++)
    //   spec->J[idl,idz] = &spectrum.J[idz][idl];
  }
  spec->J = spectrum.J;
  spec->nlw = Nlw;
  // memcpy(spec->J, spectrum.J, (Nlw+1)*atmos.Nspace * sizeof(double));
}