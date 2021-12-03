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

#define COMMENT_CHAR        "#"

#define RAY_INPUT_FILE      "ray.input"
#define ASCII_SPECTRUM_FILE "spectrum_%4.2f.asc"


/* --- Function prototypes --                          -------------- */

void get_spec(char *argv[], double muz, int pID, int idz, double perturbation);
// void get_spec(char *argv[], double muz);
void perturb_parameter(int idp, int idz, double perturbation);

/* --- Global variables --                             -------------- */

enum Topology topology = ONE_D_PLANE;

Atmosphere atmos;
Geometry geometry;
Spectrum spectrum;
ProgramStats stats;
InputData input;
CommandLine commandline;
char messageStr[MAX_LINE_SIZE];

/* ------- begin -------------------------- solveray.c -------------- */

int main(int argc, char *argv[])
{
  register int n, k, la;

  char    rayFileName[14], inputLine[MAX_LINE_SIZE], ascFilename[18];
  bool_t  result, exit_on_EOF, to_obs, initialize, crosscoupling,
          analyze_output, equilibria_only;
  int     Nspect, Nread, Nrequired, checkPoint, *wave_index = NULL;
  double  muz, *S, *chi, *J;
  FILE   *fp_out, *fp_ray, *fp_stokes, *fp_out_asc;
  XDR     xdrs;
  ActiveSet *as;

  setOptions(argc, argv);
  // getCPU(0, TIME_START, NULL);
  SetFPEtraps();

  /* --- Read input data and initialize --             -------------- */

  readInput();
  spectrum.updateJ = FALSE;

  /* --- Read direction cosine for ray --              -------------- */

  if ((fp_ray = fopen(RAY_INPUT_FILE, "r")) == NULL) {
    sprintf(messageStr, "Unable to open inputfile %s", RAY_INPUT_FILE);
    Error(ERROR_LEVEL_2, argv[0], messageStr);
  }
  
  getLine(fp_ray, COMMENT_CHAR, inputLine, exit_on_EOF=TRUE);
  Nread = sscanf(inputLine, "%lf", &muz);
  checkNread(Nread, Nrequired=1, argv[0], checkPoint=1);

  if (muz <= 0.0  ||  muz > 1.0) {
    sprintf(messageStr,
      "Value of muz = %f does not lie in interval <0.0, 1.0]\n", muz);
    Error(ERROR_LEVEL_2, argv[0], messageStr);
  }

  if (input.StokesMode == FIELD_FREE ||
      input.StokesMode == POLARIZATION_FREE) {
    input.StokesMode = FULL_STOKES;
  }

  get_spec(argv, muz, 0, -1, 0);

  int nz=atmos.Nspace, nw=spectrum.Nspect;
  int idz, idw, ids, idp;
  double spec_plus[nw][4], spec_minus[nw][4];
  long double RF[nz][nw][4][6];
  // perturbations = 1K, 1m/s, 1m/s, 1G, 1e-3rad, 1e-3rad
  double perturbation[6] = {1.0, 1.0, 1.0, 1e-4, 0.01, 0.01};
  FILE *fp_rf;

  int rf_par_flag[6] = {0,0,0,0,0,0};

  if (input.rf_temp) rf_par_flag[0] = 1;
  if (input.rf_vz) rf_par_flag[1] = 2;
  if (input.rf_vmic) rf_par_flag[2] = 3;
  if (input.rf_mag) rf_par_flag[3] = 4;
  if (input.rf_gamma) rf_par_flag[4] = 5;
  if (input.rf_chi) rf_par_flag[5] = 6;

  for (idp=0; idp<6; idp++) {
    if (rf_par_flag[idp]>0) {
      for (idz=0; idz<nz; idz++) {
        
        // positive perturbation
        get_spec(argv, muz, rf_par_flag[idp], idz, perturbation[idp]);
          
        for (idw=0; idw<nw; idw++) {
          spec_plus[idw][0] = spectrum.I[0][idw];
          spec_plus[idw][1] = spectrum.Stokes_Q[0][idw];
          spec_plus[idw][2] = spectrum.Stokes_U[0][idw];
          spec_plus[idw][3] = spectrum.Stokes_V[0][idw];
        }

        // negative perturbation
        // if we have angles, we do forward differences
        // else we do central differences
        if ((idp==4) || (idp==5)) {
          get_spec(argv, muz, rf_par_flag[idp], idz, 0);
        }
        else {
          get_spec(argv, muz, rf_par_flag[idp], idz, -perturbation[idp]);
        }
        
        for (idw=0; idw<nw; idw++) {
          spec_minus[idw][0] = spectrum.I[0][idw];
          spec_minus[idw][1] = spectrum.Stokes_Q[0][idw];
          spec_minus[idw][2] = spectrum.Stokes_U[0][idw];
          spec_minus[idw][3] = spectrum.Stokes_V[0][idw];
        }

        // return perturbation
        // perturb_parameter(rf_par_flag[idp], idz, perturbation[idp]);
        
        if ((idp==4) || (idp==5)) {
          for (idw=0; idw<nw; idw++) {
            for (ids=0; ids<4; ids++) {
              RF[idz][idw][ids][idp] = (spec_plus[idw][ids] - spec_minus[idw][ids])/perturbation[idp];
            }
            // printf("RF = %15.13Le\n", RF[idz][idw][0][idp]);
          }
        }
        else {
          for (idw=0; idw<nw; idw++) {
            for (ids=0; ids<4; ids++) {
              RF[idz][idw][ids][idp] = (spec_plus[idw][ids] - spec_minus[idw][ids])/2/perturbation[idp];
            }
          }
        }

      } // nz for loop end
    } // rf_par_flag if end
    else {
      for (idz=0; idz<nz; idz++) {
        for (idw=0; idw<nw; idw++) {
          for (ids=0; ids<4; ids++) {
            RF[idz][idw][ids][idp] = -99;
          }
        }
      }
    } // end else
  } // idp for loop end

  fp_rf = fopen(input.rfs_output, "w");
  
  for (idp=0; idp<6; idp++) {
    for (idz=0; idz<nz; idz++) {
      for (idw=0; idw<nw; idw++) {
        for (ids=0; ids<4; ids++) {
          fprintf(fp_rf, "%15.13Le ", RF[idz][idw][ids][idp]);
        }
      }
    }
    fprintf(fp_rf, "\n");
  }
  fclose(fp_rf);

  return 0;
}

void perturb_parameter(int idp, int idz, double perturbation)
{
  switch (idp)
  {
    case 1:
      // printf("Temperature\n");
      atmos.T[idz] += perturbation;
      break;
    case 2:
      // printf("Vertical velocity\n");
      geometry.vel[idz] += perturbation;
      break;
    case 3:
      // printf("Turbulent velocity\n");
      atmos.vturb[idz] += perturbation;
      break;
    case 4:
      // printf("Magnetic field\n");
      atmos.B[idz] += perturbation;
      break;
    case 5:
      // printf("Inclination\n");
      atmos.gamma_B[idz] += perturbation;
      break;
    case 6:
      // printf("Azimuth\n");
      atmos.chi_B[idz] += perturbation;
    default:
      break;
  }
}

void get_spec(char *argv[], double muz, int pID, int idz, double perturbation)
{
  bool_t analyze_output, equilibria_only;

  /* --- Read input data for atmosphere --             -------------- */

  // getCPU(1, TIME_START, NULL);
  MULTIatmos(&atmos, &geometry);

  perturb_parameter(pID, idz, perturbation);
 
  /* --- redefine geometry for just this one ray --    -------------- */
  
  atmos.Nrays = geometry.Nrays = 1;
  geometry.muz[0] = muz;
  geometry.mux[0] = sqrt(1.0 - SQ(geometry.muz[0]));
  geometry.muy[0] = 0.0;
  geometry.wmu[0] = 1.0;
  if (atmos.Stokes) Bproject();

  input.startJ = OLD_J;

  readAtomicModels();
  readMolecularModels();
  SortLambda();

  getBoundary(&geometry);

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
    convertScales(&atmos, &geometry);
    
    getProfiles();
    initSolution();
    initScatter();

    // getCPU(1, TIME_POLL, "Total initialize");

    /* --- Solve radiative transfer equations --         -------------- */

    solveSpectrum(FALSE, FALSE);
}

/* ------- end ---------------------------- solveray.c -------------- */
