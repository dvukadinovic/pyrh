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
#include "../background.h"
#include "../statistics.h"
#include "../inputs.h"
#include "../error.h"
#include "../../headers/xdr.h"
#include "../constant.h"

#include "pyrh_compute1dray.h"

/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

// enum Topology topology = ONE_D_PLANE;

// extern Atmosphere atmos;
// extern Geometry geometry;
// extern InputData input;
// extern Spectrum spectrum;


/* --- Global variables --- */

extern Atmosphere atmos;
extern Geometry geometry;
extern Spectrum spectrum;
extern ProgramStats stats;
extern InputData input;
extern CommandLine commandline;
extern char messageStr[MAX_LINE_SIZE];

// functions declaration
int _getnumber(int* z);
void _solveray(double muz, mySpectrum *spec);

int _getnumber(int* z)
{
  int a = 1;
  printf("z = %d\n", z[0]);
  return a+z[0];
}

/* ------- begin -------------------------- solveray.c -------------- */

void _solveray(double muz, mySpectrum *spec)
{
  bool_t  result, exit_on_EOF, to_obs, initialize, crosscoupling,
          analyze_output, equilibria_only;
  
  /* --- Solve radiative transfer equations --         -------------- */

  if (input.solve_NLTE){
    input.startJ = OLD_J;
    spectrum.updateJ = FALSE;
    input.limit_memory = FALSE;

    atmos.Nrays = geometry.Nrays = 1;
    geometry.muz[0] = muz;
    geometry.mux[0] = sqrt(1.0 - SQ(geometry.muz[0]));
    geometry.muy[0] = 0.0;
    geometry.wmu[0] = 1.0;
    if (atmos.Stokes) Bproject();

    if (atmos.moving || input.StokesMode) {
      Background(analyze_output=FALSE, equilibria_only=FALSE);
    } else {
      Background(analyze_output=FALSE, equilibria_only=TRUE);
    }

    getProfiles();
    // spectrum.J is already filled with correct values; do we need this initSolution() here? 
    // It seems not (tested on H in ACTIVE state; everything was exactly the same as from original RH)
    // initSolution(FALSE);
    initScatter();
    
    solveSpectrum(FALSE, FALSE);
  }
  
  spec->sQ = NULL;
  spec->sU = NULL;
  spec->sV = NULL;
  spec->J = NULL;
  spec->J20 = NULL;

  int Nlw = spec->nlw;
  Nlw -= 1;

  spec->stokes = 0;
  spec->lam = (double *) malloc(Nlw * sizeof(double));
  spec->sI = (double *) malloc(Nlw * sizeof(double));
  spec->sQ = (double *) malloc(Nlw * sizeof(double));
  spec->sU = (double *) malloc(Nlw * sizeof(double));
  spec->sV = (double *) malloc(Nlw * sizeof(double));
  if (input.get_atomic_rfs) spec->rfs = matrix_double(Nlw, input.n_atomic_pars);
  // spec->J  = matrix_double(Nlw+1, atmos.Nspace);

  int index=0;
  double tmp;
  
  for (int idl=0; idl<Nlw+1; idl++){
    // skip referent wavelength
    if (spectrum.lambda[idl]==atmos.lambda_ref) continue;
    // vacuum_to_air(1, &spectrum.lambda[idl], &tmp);
    // spec->lam[index] = tmp;
    spec->lam[index] = spectrum.lambda[idl];
    spec->sI[index] = spectrum.I[idl][0];
    if (atmos.Stokes){
      spec->sQ[index] = spectrum.Stokes_Q[idl][0];
      spec->sU[index] = spectrum.Stokes_U[idl][0];
      spec->sV[index] = spectrum.Stokes_V[idl][0];
      spec->stokes = 1;
    }
    // free_as(idl, FALSE);
    if (input.get_atomic_rfs){
      for (int idp=0; idp<input.n_atomic_pars; idp++){
        spec->rfs[index][idp] = atmos.atomic_rfs[idl][0][idp];
      }
    }
    index += 1;
  }
  // for (int idl=0; idl<spectrum.Nspect-3; idl++){
  //   printf(" -- %d :: ", idl);
  //   // free_as(idl, FALSE);
  // }

  spec->nlw = Nlw;

  // deallocate Stokes spectrum
  if (spectrum.lambda!=NULL) free(spectrum.lambda);
  if (spectrum.I!=NULL) freeMatrix((void **) spectrum.I);
  if (spectrum.Stokes_Q!=NULL) freeMatrix((void **) spectrum.Stokes_Q);
  if (spectrum.Stokes_U!=NULL) freeMatrix((void **) spectrum.Stokes_U);
  if (spectrum.Stokes_V!=NULL) freeMatrix((void **) spectrum.Stokes_V);

  // deallocate J and J20
  if (spectrum.J!=NULL) freeMatrix((void **) spectrum.J);
  if (input.backgr_pol){
    if (spectrum.J20!=NULL) freeMatrix((void **) spectrum.J20);
  }

  Atom *atom;

  if (input.get_populations){
    spec->atom_pops = malloc(atmos.Nactiveatom * sizeof(AtomPops));
    spec->Nactive_atoms = atmos.Nactiveatom;
    for (int nact = 0;  nact < atmos.Nactiveatom;  nact++) {
      atom = atmos.activeatoms[nact];

      // spec->atom_pops[nact].ID = atom->ID;
      strcpy(&spec->atom_pops[nact].ID[0], &atom->ID[0]);
      spec->atom_pops[nact].Nlevel = atom->Nlevel;
      spec->atom_pops[nact].Nz = atmos.Nspace;
      spec->atom_pops[nact].n = atom->n;
      spec->atom_pops[nact].nstar = atom->nstar;
    } 
  }
}