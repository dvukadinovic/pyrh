/* ------- file: -------------------------- background.c ------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Wed Jul 24 12:52:46 2013 --

       --------------------------                      ----------RH-- */

/* Driving subroutine for background opacity sources.

 * Included at the moment:

  ++ Thomson scattering by free electrons

  ++ Hydrogen:
    -- Bound-free absorption and emission by neutral Hydrogen
    -- Free-free absorption and emission by neutral Hydrogen
    -- Rayleigh scattering by neutral Hydrogen and Helium
    -- Rayleigh scattering by molecular Hydrogen (H2)
    -- Bound-free absorption and emission by Hminus (H^-)
    -- Free-free absorption and emission by Hminus (H + e)
    -- Free-free absorption and emission by H2minus (H2 + e)
    -- Free-free absorption and emission by H2plus (H + H^+)

  ++ Metals:
    -- Bound-free absorption and emission from metals specified in
       file background.input.
    -- Bound-bound absorption and emission from metals specified in
       file background.input.
    -- LTE Bound-bound absorption and emission all elements from
       a Kurucz line list.

  ++ Molecules:
    -- Chemical equilibrium is calculated for molecules specified
       in file background.input and populations of constituent atoms
       are reduced accordingly.
    -- molecular opacities (LTE) may be taken into account by specifying
       data files with transition lists in the molecular input files.
    -- Bound-free absorption and emission by OH and CH molecules.

 * Atomic models are specified in atoms.input, molecules in 
   molecules.input

 * Entries for the atoms.input and molecules.input files should have
   the form, respectively:

    -------------------------------------------------------------------
   |                                                                  |
   |   Nmetal                                                         |
   |                                                                  |
   |   model file ACTIVE/PASSIVE  INITIAL_SOLUTION   population file  |
   |                          .                                       |
   |                          .                                       |
   |                                                                  |
   |   Nmolecule                                                      |
   |                                                                  |
   |   model file ACTIVE/PASSIVE  INITIAL_SOLUTION                    |
   |                 .                                                |
   |                 .                                                |
    -------------------------------------------------------------------


   Nmetal and Nmolecule are the number of metal and molecule entries.
   metalID is the two-character atomID, the next entry is either
   set to LTE or NLTE, model_file is the input file containing atomic
   data for this metal (generic atomic input data format), and
   population_file is the input file containing the NLTE population
   numbers from a previous calculation. This last entry is only read when
   the second entry is set to NLTE. 

   -- Units:
      Wavelengths are given in nm, densities in m^-3, opacities in m^2,
      and emissivities in J s^-1 Hz^-1 sr^-1.

 Note: The model atom file for hydrogen is specified in keyword.input.
       If H_LTE = TRUE is specified there LTE hydrogen populations are
       used. See: distribute_nH in the file hydrogen.c

 Note: Scattering opacity is added to total opacity after all
       contributions have been computed.

 Note: The quantities chi_ai and eta_ai store the angle-inpendent
       opacities and emissivities in case atmos.moving == TRUE.
       In static atmospheres these quantities are just mapped to
       atmos.chi_c and atmos.eta_c to save memory space.

 Note: FALSE the auxiliary output files for
       the Background Record Structure (BRS), metals, and molecules
       are NOT written. This option is used when Background is called
       from programs like solveray (formal solution along one specific
       ray) in cases with moving atmospheres (angle-dependent opacity).

 Note: If equilibria_only is set to TRUE only the electron density,
       LTE populations and collisions, and chemical equilibria are
       evaluated.

 Note: Record numbers stored in atmos.backgrrecno refer to records
       of the size atmos.Nspace. If a polarized line is present 9
       (4 + 4 + 1, no magneto-optical effects), or 12 (7 + 4 +1, with
       magneto-optical effects) records are used, otherwise 3 (1 + 1 + 1)
       records.
       If the atmosphere is moving, or if a polarized line is present
       data is stored for each angle and wavelength, otherwise data
       is stored once for each wavelength only.
       --                                              -------------- */

#include <fcntl.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "constant.h"
#include "background.h"
#include "error.h"
#include "statistics.h"
#include "inputs.h"
#include "xdr.h"

#include "pyrh_background.h"

#define COMMENT_CHAR  "#"

/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Spectrum spectrum;
extern InputData input;
extern char messageStr[];

void pyrh_Background(bool_t equilibria_only, double* total_opacity)
{
  const char routineName[] = "pyrh_Background";
  register int k, nspect, n, mu, to_obs;

  static int ne_iter = 0;
  char    inputLine[MAX_LINE_SIZE];
  bool_t  exit_on_EOF, do_fudge, fromscratch, debug = FALSE;
  int     backgrrecno, index, Nfudge, NrecStokes, recordsize;
  double *chi, *eta, *scatt, wavelength, *thomson, *chi_ai, *eta_ai, *sca_ai,
          Hmin_fudge, scatt_fudge, metal_fudge, *lambda_fudge, **fudge,
         *Bnu, *chi_c, *eta_c, *sca_c, *chip, *chip_c;
  Atom   *He;
  FILE   *fp_fudge;
  flags   backgrflags;

  getCPU(2, TIME_START, NULL);

  /* --- Set equilibrium state -------------------------------------- */

  // printf("Getting equilibrium state\n");
  // fromscratch = TRUE; // because SOLVE_NE == ONCE always for HSE!
  // if (debug) printf("Solve_ne()\n");
  // _Solve_ne(fromscratch);
  // if (debug) printf("Set LTE\n");
  // SetLTEQuantities();
  
  // printf("Setting ChemicalEquilibrium\n");
  if (input.NonICE)
    readMolecules(MOLECULAR_CONCENTRATION_FILE);
  else
    if (debug) printf("ChemEq\n");
    ChemicalEquilibrium(N_MAX_CHEM_ITER, CHEM_ITER_LIMIT);

  if (equilibria_only) {

    /* --- If we only need ne, LTE populations and collisions, and
           chemical equilibrium leave here --          -------------- */

    getCPU(2, TIME_POLL, "Total Background");
    return;
  }
    
  getCPU(3, TIME_START, NULL);

  do_fudge = input.do_fudge;
  if (do_fudge){
    Nfudge = atmos.fudge_num;
    lambda_fudge = atmos.fudge_lam;
    fudge = atmos.fudge;
  }

  /* --- Allocate temporary storage space. The quantities are used
         for the following purposes:

       - chi, eta, scatt: Get contributions to opacity, emissivity,
         and scattering opacity, respectively, from a specific process
         for a given wavelength and possibly angle.

       - chi_c, eta_c, sca_c: Collect the total opacity, emissivity
         and scattering opacity for a given wavelength and possibly
         angle.

       - chi_ai, eta_ai: Collect the angle-independent part of
         opacity and emissivity for each wavelength so that these
         need not be recalculated in an angle-dependent case.
         When the atmosphere is not moving and has no magnetic fields
         these just point to the total quantities chi_c and eta_c.

   Note: In case of magnetic fields in the atmosphere chi, eta and 
         chip, and chi_c, eta_c and chip_c contain all four Stokes
         parameters, and should be allocated a 4 and 3 times larger
         storage space, respectively.
         --                                            -------------- */

  if (atmos.Stokes)
    NrecStokes = 4;
  else
    NrecStokes = 1;

  // we are only doing it layer by layer
  int Nspace = 1;

  chi_ai = (double *) malloc(Nspace * sizeof(double));
  eta_ai = (double *) malloc(Nspace * sizeof(double));
  sca_ai = (double *) malloc(Nspace * sizeof(double));

  chi   = (double *) malloc(NrecStokes*Nspace * sizeof(double));
  eta   = (double *) malloc(NrecStokes*Nspace * sizeof(double));
  scatt = (double *) malloc(Nspace * sizeof(double));

  Bnu = (double *) malloc(Nspace * sizeof(double));

  /* --- Thomson scattering by free electrons is wavelength independent
         in non-relativistic limit so we compute it only once -- ---- */

  if (debug) printf("Get Thomson\n");
  thomson = (double *) malloc(Nspace * sizeof(double));
  Thomson(thomson);

  /* --- Check whether an atomic model is present for He -- --------- */

  He = (atmos.elements[1].model) ? atmos.elements[1].model : NULL;

  /* --- Read background files from Kurucz data file -- ------------- */

  // We do not read Kurucz lines (we used pyrh_Background() to set 
  // atmosphere in HSE only; no need for lines here).
  // if (atmos.Nrlk == 0){
  //   readKuruczLines(input.KuruczData);
  // }
  
  // if (atmos.Nrlk > 0) {
  //   qsort(atmos.rlk_lines, atmos.Nrlk, sizeof(RLK_Line), rlk_ascend);
  // }
  /* --- Allocate memory for the boolean array that stores whether
         a wavelength overlaps with a Bound-Bound transition in the
         background, or whether it is polarized --     -------------- */

  atmos.backgrflags = (flags *) malloc(spectrum.Nspect * sizeof(flags));
  for (nspect = 0;  nspect < spectrum.Nspect;  nspect++) {
    atmos.backgrflags[nspect].hasline = FALSE;
    atmos.backgrflags[nspect].ispolarized = FALSE;
  }
  /* --- Allocate memory for the list of record numbers that specifies
         for each wavelength where to find the background opacity,
         scattering opacity, and emissivity --         -------------- */

  // backgrrecno = 0;

  // if (atmos.moving || atmos.Stokes) {
  //   atmos.backgrrecno = 
  //     (int *) malloc(2*spectrum.Nspect*atmos.Nrays * sizeof(int));
  // } else
  //   atmos.backgrrecno = (int *) malloc(spectrum.Nspect * sizeof(int));

  /* --- Go through the spectrum and add the different opacity and
         emissivity contributions. This is the main loop --  -------- */

  k = 0;
  int layer = atmos.active_layer;

  if (debug) printf("Entering the lambda loop\n");

  for (nspect = 0;  nspect < spectrum.Nspect;  nspect++) {
    wavelength = spectrum.lambda[nspect];

    /* --- The Planck function at this wavelength --   -------------- */

    if (debug) printf("Get Planck function\n");
    Planck(atmos.Nspace, atmos.T, wavelength, Bnu, layer);
    
    /* --- Initialize the flags for this wavelength -- -------------- */

    atmos.backgrflags[nspect].hasline     = FALSE;
    atmos.backgrflags[nspect].ispolarized = FALSE;

    /* --- Initialize angle-independent quantities --  -------------- */

    chi_ai[k] = 0.0;
    eta_ai[k] = 0.0;
    sca_ai[k] = thomson[k];
    total_opacity[layer] = 0.0;

    /* --- Negative hydrogen ion, bound-free and free-free -- ------- */

    if (Hminus_bf(wavelength, chi, eta)) {
      if (debug) printf("Get Hminus_bf opacity\n");
    	chi_ai[k] += chi[k];
    	eta_ai[k] += eta[k];
      total_opacity[layer] += chi[k];
    }

    if (Hminus_ff(wavelength, chi)) {
      chi_ai[k] += chi[k];
      eta_ai[k] += chi[k] * Bnu[k];
      total_opacity[layer] += chi[k];
    }

    /* --- Opacity fudge factors, applied to Hminus opacity -- ------ */

    if (do_fudge) {
      Linear(Nfudge, lambda_fudge, fudge[0],
	     1, &wavelength, &Hmin_fudge, FALSE);
      chi_ai[k] *= Hmin_fudge;
    	eta_ai[k] *= Hmin_fudge;
      total_opacity[layer] *= Hmin_fudge;
    }
  
    /* --- Opacities from bound-free transitions in OH and CH -- ---- */

    // if (OH_bf_opac(wavelength, chi, eta)) {
    //   if (debug) printf("Get OH_bf opacity\n");
    // 	chi_ai[k] += chi[k];
    // 	eta_ai[k] += eta[k];
    //   total_opacity[layer] += chi[k];
    // }
    
    // if (CH_bf_opac(wavelength, chi, eta)) {
    //   if (debug) printf("Get CH_bf opacity\n");
    // 	chi_ai[k] += chi[k];
    // 	eta_ai[k] += eta[k];
    //   total_opacity[layer] += chi[k];
    // }
    
     /* --- Neutral Hydrogen Bound-Free and Free-Free --  ------------ */

    if (Hydrogen_bf(wavelength, chi, eta)) {
      if (debug) printf("Get Hydrogen_bf opacity\n");
    	chi_ai[k] += chi[k];
    	eta_ai[k] += eta[k];
      total_opacity[layer] += chi[k];
    }

    if (debug) printf("Get Hydrogen_ff opacity\n");
    Hydrogen_ff(wavelength, chi);
    chi_ai[k] += chi[k];
    eta_ai[k] += chi[k] * Bnu[k];
    total_opacity[layer] += chi[k];
    
    /* --- Rayleigh scattering by neutral hydrogen --  -------------- */

    // if (Rayleigh(wavelength, atmos.H, scatt)) {
    //   if (debug) printf("Get Rayleigh opacity\n");
	   //  sca_ai[k]  += scatt[k];
    // }

    /* --- Rayleigh scattering by neutral helium --    -------------- */
    // if (He && Rayleigh(wavelength, He, scatt)) {
    //   if (debug) printf("Get Rayleigh opacity\n");
	   //  sca_ai[k]  += scatt[k];
    // }
    /* --- Absorption by H + H^+ (referred to as H2plus free-free) -- */

    if (H2plus_ff(wavelength, chi)) {
      if (debug) printf("Get H2plus_ff opacity\n");
    	chi_ai[k] += chi[k];
    	eta_ai[k] += chi[k] * Bnu[k];
      total_opacity[layer] += chi[k];
    }

    /* --- Rayleigh scattering and free-free absorption by
           molecular hydrogen --                       -------------- */

    // if (Rayleigh_H2(wavelength, scatt)) {
    //   if (debug) printf("Get Rayleigh_H2 opacity\n");
	   //  sca_ai[k]  += scatt[k];
    // }
    
    // if (H2minus_ff(wavelength, chi)) {
    //   if (debug) printf("Get H2minus_ff opacity\n");
    // 	chi_ai[k] += chi[k];
    // 	eta_ai[k] += chi[k] * Bnu[k];
    //   total_opacity[layer] += chi[k];
    // }
    
    /* --- Bound-Free opacities due to ``metals'' --   -------------- */

    // metal_fudge = 1.0;
    // if (do_fudge) {
    //   Linear(Nfudge, lambda_fudge, fudge[2],
	   //   1, &wavelength, &metal_fudge, FALSE);
    // }
    
    /* --- Note: Hydrogen bound-free opacities are calculated in
           routine Hydrogen_bf --                      -------------- */

    // printf("layer = %d\n", atmos.active_layer);
    // Metal_bf(wavelength, atmos.Natom-1, atmos.atoms+1, chi, eta);
    // chi_ai[k] += chi[k] * metal_fudge;
    // eta_ai[k] += eta[k] * metal_fudge;

    /* --- Add the scattering opacity to the absorption part to store
           the total opacity --                        -------------- */

    scatt_fudge = 1.0;
    if (do_fudge) {
      Linear(Nfudge, lambda_fudge, fudge[1],
	     1, &wavelength, &scatt_fudge, FALSE);
    }
    sca_ai[k] *= scatt_fudge;
    chi_ai[k] += sca_ai[k];
    total_opacity[layer] += sca_ai[k];
  }

  getCPU(3, TIME_POLL, "Background Opacity");

  /* --- Free the temporary space allocated in the ff routines -- --- */

  // Hminus_ff(0.0, NULL);
  // printf("Cleaning!\n");
  // H2minus_ff(0.0, NULL);
  // printf("Cleaning!\n");
  // H2plus_ff(0.0, NULL);

  // printf("Cleaning!\n");
  free(chi);
  // printf("Cleaning!\n");
  free(eta);
  // printf("Cleaning!\n");  
  free(scatt);
  // printf("Cleaning!\n");
  free(Bnu);
  // printf("Cleaning!\n");
  free(thomson);

  // printf("Cleaning!\n");
  free(chi_ai);
  free(eta_ai);
  free(sca_ai);
  
  getCPU(2, TIME_POLL, "Total Background");
}

void get_ne(bool_t fromscratch){
  const char routineName[] = "get_ne";
  register int k, n, j;

  int     Nmaxstage, niter;
  double *fjk, *dfjk, error, ne_old, akj, sum, PhiH, C1, Uk,
    dne, dnemax, *np, PhiHmin;
  Element tmp;

  tmp = atmos.elements[0];

  getCPU(3, TIME_START, NULL);

  C1 = (HPLANCK/(2.0*PI*M_ELECTRON)) * (HPLANCK/KBOLTZMANN);

  /* --- Figure out the largest array size needed so that we do not
         have to allocate and free memory all the time -- ----------- */

  Nmaxstage = 0;
  for (n = 0;  n < atmos.Nelem;  n++)
    Nmaxstage = MAX(Nmaxstage, atmos.elements[n].Nstage);
  fjk  = (double *) malloc(Nmaxstage * sizeof(double));
  dfjk = (double *) malloc(Nmaxstage * sizeof(double));

  int layer = atmos.active_layer;

  np = atmos.H->n[atmos.H->Nlevel-1];
  if (fromscratch) {
    /* --- Get the initial solution from ionization of H only -- -- */
    if (atmos.H_LTE) {
      // Uk = getKuruczpf(&atmos.elements[0], 0, layer);
      // Forced Uk for Hydrogen to 0 because on MPS clusters this has 
      // value of 45 for some reason... 
      // Partition functions are the same, temperature is the same, nHtot is
      // the same. Not still clear why is there a difference. Synthesis works fine,
      // only problem pertains in HSE computation for some reason.
      Uk = 0;
      PhiH = 0.5 * pow(C1/atmos.T[layer], 1.5) *
         exp(Uk + atmos.elements[0].ionpot[0]/(KBOLTZMANN*atmos.T[layer]));
      ne_old = (sqrt(1.0 + 4.0*atmos.nHtot[layer]*PhiH) - 1.0) / (2.0*PhiH);
      //if (layer==0) printf("PhiH = %e | T = %f | Uk = %e | Ej = %e\n", PhiH, atmos.T[layer], Uk, atmos.elements[0].ionpot[0]);
    } else
       ne_old = np[layer];
       /* --- Copy into ne as well to calculate first fij and dfij - - */
       atmos.ne[layer] = ne_old;
  } else {
    /* --- Use original electron density as starting guess -- ----- */
    ne_old = atmos.ne[layer];
  }

  niter = 0;
  while (niter < N_MAX_ELECTRON_ITERATIONS) {
    error = ne_old / atmos.nHtot[layer];
    sum   = 0.0;
    for (n = 0;  n < atmos.Nelem;  n++) {
      getfjk(&atmos.elements[n], ne_old, layer, fjk, dfjk);
      /* --- Contribution from Hminus --             -------------- */
      if (n == 0) {
        PhiHmin = 0.25*pow(C1/atmos.T[layer], 1.5) *
            exp(E_ION_HMIN / (KBOLTZMANN * atmos.T[layer]));
        error += ne_old * fjk[0] * PhiHmin;
        sum   -= (fjk[0] + ne_old * dfjk[0]) * PhiHmin;
      }
      for (j = 1;  j < atmos.elements[n].Nstage;  j++) {
        akj = atmos.elements[n].abund * j;
        error -= akj * fjk[j];
        sum   += akj * dfjk[j];
      }
    }

    atmos.ne[layer] = ne_old -
        atmos.nHtot[layer] * error / (1.0 - atmos.nHtot[layer] * sum);
    dne = fabs((atmos.ne[layer] - ne_old)/ne_old);
    ne_old = atmos.ne[layer];
  
    if (dne <= MAX_ELECTRON_ERROR) break;
    niter++;
  }

  if (dne > MAX_ELECTRON_ERROR) {
    sprintf(messageStr, "Electron density iteration not converged:\n"
      " spatial location: %d, temperature: %6.1f [K], \n"
      " density: %9.3E [m^-3],\n dnemax: %9.3E\n",
      k, atmos.T[k], atmos.nHtot[k], dne);
    Error(WARNING, routineName, messageStr);
  }

  free(fjk);  free(dfjk);

  getCPU(3, TIME_POLL, "Electron density");
}

void get_Nm_total(double* Nm, int k){
  // total number of molecular density in atmosphere at k-th depth
  Nm[k] = 0;
  for (int n=0;  n<atmos.Nmolecule; n++){
    Nm[k] += atmos.molecules[n].n[k];
  }
}
