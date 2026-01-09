/* ------- file: -------------------------- stokesopac.c ------------

       Version:       rh1.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Mon Feb 23 15:51:02 2004 --

       --------------------------                      ----------RH-- */

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "inputs.h"

/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Spectrum spectrum;
extern InputData input;
extern char messageStr[];


/* ------- begin -------------------------- StokesK.c --------------- */

void StokesK(int nspect, int k, double chi_I, double K[4][4])
{
  register int i, j;

  ActiveSet *as;

  /* --- Return the elements of the reduced 4x4 Stokes propagation
         matrix K', which is defined as:

           =           =         =              =
           K' = (chi_c*1 + chi_l*Phi) / chi_I - 1,

	   for wavelength# nspect, spatial point k, and ray mu.

    See: Rees, Murphy, & Durrant, 1989, ApJ 339, 1093-1106.

         --                                            -------------- */
  as = &spectrum.as[nspect];

  for (j = 0;  j < 4;  j++)
    for (i = 0;  i < 4;  i++) K[j][i] = 0.0;

  if (containsPolarized(as)) {
    K[0][1] = as->chi[atmos.Nspace + k];
    K[0][2] = as->chi[2*atmos.Nspace + k];
    K[0][3] = as->chi[3*atmos.Nspace + k];

    if (input.magneto_optical) {
      K[1][2] = as->chip[2*atmos.Nspace + k];
      K[1][3] = as->chip[atmos.Nspace + k];
      K[2][3] = as->chip[k];
    }
  }
  if (atmos.backgrflags[nspect].ispolarized) {
    K[0][1] += as->chi_c[atmos.Nspace + k];
    K[0][2] += as->chi_c[2*atmos.Nspace + k];
    K[0][3] += as->chi_c[3*atmos.Nspace + k];

    if (input.magneto_optical) {
      K[1][2] += as->chip_c[2*atmos.Nspace + k];
      K[1][3] += as->chip_c[atmos.Nspace + k];
      K[2][3] += as->chip_c[k];
    }
  }
  /* --- Divide by Stokes I opacity and fill lower diagonal part -- - */

  for (j = 0;  j < 3;  j++) {
    for (i = j+1;  i < 4;  i++) {
      K[j][i] /= chi_I;
      K[i][j]  = K[j][i];
    }
  }
  /* --- Anti-symmetric magneto-optical elements --    -------------- */

  if (input.magneto_optical) {
    K[1][3] *= -1.0;
    K[2][1] *= -1.0;
    K[3][2] *= -1.0;
  }
}
/* ------- end ---------------------------- StokesK.c --------------- */

void Gvec(int nspect, int k, double chi_I, double **I, double **G)
{
  register int i, j, idp;
  double mat[4][4];

  // reset everything to zero
  for (idp=0; idp<input.n_atomic_pars; idp++) {
    for (j = 0;  j < 4;  j++){
      G[idp][j] = 0.0;
    }
  }

  for (idp=0; idp<input.n_atomic_pars; idp++){
    // memset(&mat[0][0], 0, 16*sizeof(double));
    for (j=0; j<4; j++)
      for (i=0; i<4; i++)
        mat[j][i] = 0.0;

    // get the opacity perturbation matrix
    for (j = 0; j<4; j++) mat[j][j] = spectrum.dchi_c_lam[nspect][k][idp];
    
    if (atmos.backgrflags[nspect].ispolarized) {
      mat[0][1] = spectrum.dchi_Q[nspect][k][idp];
      mat[0][2] = spectrum.dchi_U[nspect][k][idp];
      mat[0][3] = spectrum.dchi_V[nspect][k][idp];
      
      if (input.magneto_optical) {
        mat[1][2] = spectrum.dchip_V[nspect][k][idp];
        mat[1][3] = spectrum.dchip_U[nspect][k][idp];
        mat[2][3] = spectrum.dchip_Q[nspect][k][idp];
      }
    }

    // fill lower diagonal part
    for (j = 0;  j < 3;  j++) {
      for (i = j+1;  i < 4;  i++) {
        mat[i][j] = mat[j][i];
      }
    }
    /* --- Anti-symmetric magneto-optical elements --    -------------- */

    if (input.magneto_optical) {
      mat[1][3] *= -1.0;
      mat[2][1] *= -1.0;
      mat[3][2] *= -1.0;
    }

    // multiply opacity perturbation matrix with Stokes vector I
    for (j=0; j<4; j++){
      for (i=0; i<4; i++){
        G[idp][j] += mat[j][i] * I[i][k];
      }
      G[idp][j] *= -1;
    }

    // add emissivity perturbation
    G[idp][0] += spectrum.deta_c_lam[nspect][k][idp];

    if (atmos.backgrflags[nspect].ispolarized) {
      G[idp][1] += spectrum.deta_Q[nspect][k][idp];
      G[idp][2] += spectrum.deta_U[nspect][k][idp];
      G[idp][3] += spectrum.deta_V[nspect][k][idp];
    }

    // divide everything with the total opacity in Stokes I
    for (j=0; j<4; j++){
      G[idp][j] /= chi_I;
    }
  }

}