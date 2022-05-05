/* ------- file: -------------------------- multiatmos.c ------------

       Version:       rh2.0, 1-D plane-parallel
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Mon May 21 13:51:30 2018 --

       --------------------------                      ----------RH-- */

/* --- Reads atmospheric model in MULTI format. --     -------------- */


#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "spectrum.h"
#include "background.h"
#include "constant.h"
#include "error.h"
#include "inputs.h"
#include "statistics.h"
#include "xdr.h"


#define MULTI_COMMENT_CHAR  "*"
#define N_HYDROGEN_MULTI     6


/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern Spectrum spectrum;
extern InputData input;
extern char messageStr[];


/* ------- begin -------------------------- MULTIatmos.c ------------ */

void MULTIatmos(Atmosphere *atmos, Geometry *geometry)
{
  const char routineName[] = "MULTIatmos";
  register int k, n, mu;

  char    scaleStr[20], inputLine[MAX_LINE_SIZE], *filename;
  bool_t  exit_on_EOF, enhanced_atmos_ID = FALSE;
  int     Nread, Ndep, Nrequired, checkPoint;
  double *dscale, turbpress, turbelecpress, nbaryon, meanweight;

  getCPU(2, TIME_START, NULL);

  /* --- Get abundances of background elements --        ------------ */
 
  readAbundance(atmos);

  atmos->NHydr = N_HYDROGEN_MULTI;

  /* --- Boundary condition at TOP of atmosphere --      ------------ */

  if (strcmp(input.Itop, "none"))
    geometry->vboundary[TOP] = IRRADIATED;
  else 
    geometry->vboundary[TOP] = ZERO;

  /* --- Boundary condition at BOTTOM of atmosphere --   ------------ */

  geometry->vboundary[BOTTOM] = THERMALIZED;

  strcpy(atmos->ID, "dummyatmos");
  atmos->gravity = 4.4;

  /* --- Keep duplicates of some of the geometrical quantities in
         Atmos structure --                            -------------- */

  atmos->Ndim = 1;
  atmos->N = (int *) malloc(atmos->Ndim * sizeof(int));
  atmos->Nspace = Ndep = geometry->Ndep;
  atmos->N[0] = Ndep;

  atmos->gravity = POW10(atmos->gravity) * CM_TO_M;

  /* --- Allocate space for arrays that define structure -- --------- */

  geometry->tau_ref = (double *) malloc(Ndep * sizeof(double));
  geometry->cmass   = (double *) malloc(Ndep * sizeof(double));
  geometry->height  = (double *) malloc(Ndep * sizeof(double));

  getAngleQuad(geometry);
  atmos->wmu = geometry->wmu;

  getCPU(2, TIME_POLL, "Read Atmosphere");
}
/* ------- end ---------------------------- MULTIatmos.c ------------ */

/* ------- begin -------------------------- convertScales.c --------- */

void convertScales(Atmosphere *atmos, Geometry *geometry)
{
  register int k;

  bool_t hunt;
  int    ref_index, Ndep = geometry->Ndep;
  double *rho, *height, *cmass, *tau_ref, h_zero, unity;
  ActiveSet *as;

  /* --- Convert between different depth scales --       ------------ */

  height  = geometry->height;
  cmass   = geometry->cmass;
  tau_ref = geometry->tau_ref;

  rho = (double *) malloc(Ndep * sizeof(double));
  for (k = 0;  k < Ndep;  k++)
    rho[k] = (AMU * atmos->wght_per_H) * atmos->nHtot[k];

  /* --- Get opacity of reference wavelength --          ------------ */

  Locate(spectrum.Nspect, spectrum.lambda, atmos->lambda_ref, &ref_index);

  as = &spectrum.as[ref_index];
  alloc_as(ref_index, FALSE);
  readBackground(ref_index, 0, 0);

  /* --- Convert to missing depth scales --              ------------ */

  switch (geometry->scale) {
  case COLUMN_MASS:
    height[0] = 0.0;
    tau_ref[0] = as->chi_c[0] / rho[0] * cmass[0];
    for (k = 1;  k < Ndep;  k++) {
          height[k] = height[k-1] - 2.0*(cmass[k] - cmass[k-1]) / 
    	(rho[k-1] + rho[k]);
          tau_ref[k] = tau_ref[k-1] + 0.5*(as->chi_c[k-1] + as->chi_c[k]) *
    	(height[k-1] - height[k]);
            // height[k] = height[k-1] + 1/rho[k] * (cmass[k] - cmass[k-1]);
            // tau_ref[k] = tau_ref[k-1] + as->chi_c[k]/rho[k] * (cmass[k] - cmass[k-1]);
    }
    break;
  case TAU500:
    height[0] = 0.0;
    cmass[0]  = (tau_ref[0] / as->chi_c[0]) * rho[0];
    for (k = 1;  k < Ndep;  k++) {
          height[k] = height[k-1] - 2.0 * (tau_ref[k] - tau_ref[k-1]) / 
    	(as->chi_c[k-1] + as->chi_c[k]);
          cmass[k]  = cmass[k-1]  + 0.5*(rho[k-1] + rho[k]) *
    	(height[k-1] - height[k]);
          // height[k] = height[k-1] + 1/as->chi_c[k] * (tau_ref[k] - tau_ref[k-1]);
          // cmass[k] = cmass[k-1] + rho[k]/as->chi_c[k] * (tau_ref[k] - tau_ref[k-1]);
        // printf("%e | %f\n", as->chi_c[k], tau_ref[k]);
    }
    break;
  case GEOMETRIC:
    cmass[0] = (atmos->nHtot[0] * atmos->totalAbund + atmos->ne[0]) *
      (KBOLTZMANN * atmos->T[0] / atmos->gravity);
    tau_ref[0] = 0.5 * as->chi_c[0] * (height[0] - height[1]);
    if (tau_ref[0] > 1.0) tau_ref[0] = 0.0;
    for (k = 1;  k < Ndep;  k++) {
      cmass[k]  = cmass[k-1]  + 0.5*(rho[k-1] + rho[k]) *
	(height[k-1] - height[k]);
      tau_ref[k] = tau_ref[k-1] + 0.5*(as->chi_c[k-1] + as->chi_c[k]) *
	(height[k-1] - height[k]);
    }
    break;
  }
  free_as(ref_index, FALSE);

  if ((geometry->scale == COLUMN_MASS) || (geometry->scale == TAU500)) {
    unity = 1.0;
    Linear(Ndep, tau_ref, height, 1, &unity, &h_zero, hunt=FALSE);
    for (k = 0;  k < Ndep;  k++)
      {
        height[k] = height[k] - h_zero;
        // printf("%e\n", height[k]/1e3);
      }
  }

  free(rho);
}
/* ------- end ---------------------------- convertScales.c --------- */

/* ------- begin -------------------------- getBoundary.c ----------- */

void getBoundary(Geometry *geometry)
{
  const char routineName[] = "getBoundary";
  register int la;

  bool_t result = TRUE;
  FILE  *fp_Itop;
  XDR    xdrs;

  switch (geometry->vboundary[TOP]) {
  case ZERO: break;
  case THERMALIZED: break;
  case IRRADIATED:

    sprintf(messageStr, "\n -- reading irradiance input file: %s\n\n",
	    input.Itop);
    Error(MESSAGE, NULL, messageStr);

    geometry->Itop = matrix_double(spectrum.Nspect, geometry->Nrays);

    /* --- Open input file for irradiation at TOP --     -------------- */

    if ((fp_Itop = fopen(input.Itop, "r")) == NULL) {
      sprintf(messageStr, "Unable to open inputfile %s", input.Itop);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
    xdrstdio_create(&xdrs, fp_Itop, XDR_DECODE);

    result &= xdr_vector(&xdrs, (char *) geometry->Itop[0],
			 spectrum.Nspect * geometry->Nrays,
                         sizeof(double), (xdrproc_t) xdr_double);
    if (!result) {
      sprintf(messageStr,
	      "Unable to read irradiation data at TOP of atmosphere");
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
    xdr_destroy(&xdrs);
    fclose(fp_Itop);
    break;
  case REFLECTIVE:
    break;
  default:
    Error(ERROR_LEVEL_2, routineName,
	  "Invalid boundary condition at the TOP of atmosphere");
  }

  switch (geometry->vboundary[BOTTOM]) {
  case ZERO: break;
  case THERMALIZED: break;
  case IRRADIATED:
    geometry->Ibottom = matrix_double(spectrum.Nspect, geometry->Nrays);

    /* --- Infalling intensities at BOTTOM should be read here -- --- */

    Error(ERROR_LEVEL_1, routineName,
	  "Boundary condition IRRADIATED at BOTTOM not yet implemented");
    break;
  case REFLECTIVE:
    break;
  default:
    Error(ERROR_LEVEL_2, routineName,
	  "Invalid boundary condition at the BOTTOM of atmosphere");
  }
}
/* ------- end ---------------------------- getBoundary.c ----------- */
