/* ------- file: -------------------------- backgropac_xdr.c --------

       Version:       rh2.0
       Author:        Dusan Vukadinovic (vukadinovic@mps.mpg.de)
       Last modified: Tue Jun  8 20:00:00 2021 --

       --------------------------                      ----------RH-- */

/* --- Get opacity at reference wavelength -- - */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "background.h"
#include "error.h"
#include "inputs.h"
#include "xdr.h"

#include "geometry.h"
#include "background.h"
#include "error.h"
#include "statistics.h"

#include "constant.h"

#define COMMENT_CHAR  "#"

#define MAX_ELECTRON_ERROR         1.0E-2
#define N_MAX_ELECTRON_ITERATIONS  20

/* --- Function prototypes --                          -------------- */

void backgrOpac_(double *chi_c, double *scatt_c, int ki, double Hmin_fudge, double metal_fudge, double scatt_fudge);
void Solve_ne_(double *ne, bool_t fromscratch, int ki);
double getKuruczpf_(Element *element, int stage, int k);
void write_atmos();

/* --- Global variables --                             -------------- */

enum Topology topology = ONE_D_PLANE;

Atmosphere atmos;
Geometry geometry;
Spectrum spectrum;
InputData input;
ProgramStats stats;
CommandLine commandline;
char messageStr[MAX_LINE_SIZE];

int main(int argc, char *argv[])
{
  const char routineName[] = "getOpacity";
  double *chi_c, *scatt_c, *lambda_fudge, **fudge,
  		  Hmin_fudge, scatt_fudge, metal_fudge,
  		  Ntot, *rho, dN, break_me=1e-2,
  		  kappa_mean, Pg, Pe, Pg_old, Pe_old, dP, delta_tau;
  int     Nfudge, n, iter, k, MAX_ITER=20;

  char    inputLine[MAX_LINE_SIZE];
  bool_t  exit_on_EOF, do_fudge;
  FILE   *fp_fudge;

  double lambda = atmos.lambda_ref;

  /* --- Read input data and initialize --             -------------- */
  setOptions(argc, argv);
  SetFPEtraps();

  readInput();
  MULTIatmos(&atmos, &geometry);
  readAtomicModels();
  readMolecularModels();

  /*--- get opacity fudge factors from file ---*/
  if (strcmp(input.fudgeData, "none")) {
    do_fudge = TRUE;

    /* --- Read wavelength-dependent fudge factors to compensate for
           missing UV backround line haze --           -------------- */

    if ((fp_fudge = fopen(input.fudgeData, "r")) == NULL) {
      sprintf(messageStr, "Unable to open input file %s", input.fudgeData);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
    sprintf(messageStr,
	    "\n-Fudging background opacities with file\n  %s\n\n",
	    input.fudgeData);
    Error(MESSAGE, routineName, messageStr);

    getLine(fp_fudge, COMMENT_CHAR, inputLine, exit_on_EOF=TRUE);
    sscanf(inputLine, "%d", &Nfudge);
    lambda_fudge = (double *) malloc(Nfudge * sizeof(double));
    fudge = matrix_double(3, Nfudge);
    for (n = 0;  n < Nfudge;  n++) {
      getLine(fp_fudge, COMMENT_CHAR, inputLine, exit_on_EOF=TRUE);
      sscanf(inputLine, "%lf %lf %lf %lf", &lambda_fudge[n],
	     &fudge[0][n], &fudge[1][n], &fudge[2][n]);
    }
    for (n = 0;  n < 3*Nfudge;  n++) fudge[0][n] += 1.0;
    fclose(fp_fudge);
  } else
    do_fudge = FALSE;

  /* --- get fudge factors for reference wavelength --- */
  if (do_fudge)
  {
	  Linear(Nfudge, lambda_fudge, fudge[0],
		     1, &lambda, &Hmin_fudge, FALSE);
      Linear(Nfudge, lambda_fudge, fudge[1],
	     1, &lambda, &scatt_fudge, FALSE);
	  Linear(Nfudge, lambda_fudge, fudge[2],
	     1, &lambda, &metal_fudge, FALSE);
  } else{
  	Hmin_fudge = 1.0;
  	scatt_fudge = 1.0;
  	metal_fudge = 1.0;
  }

  chi_c   = (double *) malloc(atmos.Nspace * sizeof(double));
  scatt_c = (double *) malloc(atmos.Nspace * sizeof(double));
  rho     = (double *) malloc(atmos.Nspace * sizeof(double));

  /* --- iterate solution for HSE --- */

   // Initial guess for particple density at the top (nHtot). Assume pressure of 0.3 in CGS
  atmos.T[0] = 7485.78;
  Pg_old = 332;
  Pe = 1;

  Ntot = Pg_old/10/KBOLTZMANN/atmos.T[0];
  atmos.ne[0] = Pe/10/KBOLTZMANN/atmos.T[0];
  atmos.nHtot[0] = (Ntot - atmos.ne[0]) / atmos.totalAbund;

  printf("init HSE -- %e | %e\n", atmos.nHtot[0], atmos.ne[0]);

  for (iter=0;iter<MAX_ITER;++iter){
  	backgrOpac_(chi_c, scatt_c, k, Hmin_fudge, metal_fudge, scatt_fudge);
    rho[0] = (AMU * atmos.wght_per_H) * atmos.nHtot[0] * atmos.totalAbund;

    if (geometry.scale==COLUMN_MASS) {
      geometry.cmass[0] = geometry.tau_ref[0] * rho[0] / chi_c[0];
      Pg = geometry.cmass[0] * atmos.gravity;
      Ntot = Pg/KBOLTZMANN/atmos.T[0];
      dN = (Ntot - atmos.ne[0]) / atmos.totalAbund - atmos.nHtot[0];
    }
    if (geometry.scale==TAU500){
      Pg = (geometry.tau_ref[0] * atmos.gravity * rho[0] / chi_c[0]);
      Ntot = Pg/KBOLTZMANN/atmos.T[0];
      dN = (Ntot - atmos.ne[0]) / atmos.totalAbund - atmos.nHtot[0];
    }
    atmos.nHtot[0] += dN;
    if (fabs(dN/atmos.nHtot[0])<break_me)
      break;
  }
  printf(" %d -- %e | %e\n", iter, atmos.nHtot[0], atmos.ne[0]);
  Pg = (atmos.nHtot[0] * atmos.totalAbund + atmos.ne[0]) * KBOLTZMANN * atmos.T[0];
  printf("Pe = %e | Pg = %e | kappa = %e\n", atmos.ne[0]*KBOLTZMANN*atmos.T[0]*10, Pg*10, rho[0]/chi_c[0]/10);

  return 0;

  // iterate HSE solution over each depth
  for (k=1; k<atmos.Nspace; k++)
  // for (k=1; k<2; k++)
  {
  	// printf("%d\n", k);
  	// printf("init HSE -- %e\n", atmos.H->n[0][k]);
    if (geometry.scale==TAU500)
    {
        delta_tau = geometry.tau_ref[k] - geometry.tau_ref[k-1];
        dP = atmos.gravity * rho[k-1] / chi_c[k-1] * delta_tau;
    }
    if (geometry.scale==COLUMN_MASS)
    {
      dP = atmos.gravity * (geometry.cmass[k] - geometry.cmass[k-1]);
    }

    // Pg = (atmos.nHtot[k-1]*atmos.totalAbund + atmos.ne[k-1]) * KBOLTZMANN * atmos.T[k-1];
    // Ntot = (Pg + dP) / KBOLTZMANN / atmos.T[k] - atmos.ne[k];
    // atmos.nHtot[k] = Ntot / atmos.totalAbund;

  	// printf("init HSE %d --> %e\n", k, atmos.nHtot[k]);
    for (iter=0; iter<MAX_ITER; iter++)
    {
      backgrOpac_(chi_c, scatt_c, k, Hmin_fudge, metal_fudge, scatt_fudge); 
      rho[k] = (AMU * atmos.wght_per_H) * atmos.nHtot[k] * atmos.totalAbund;
      
      if (geometry.scale==TAU500)
      {
        kappa_mean = sqrt(chi_c[k] * chi_c[k-1] / rho[k] / rho[k-1]);
        dP = atmos.gravity/kappa_mean * delta_tau;
      }
      if (geometry.scale==COLUMN_MASS)
      {
      	dP = atmos.gravity * (geometry.cmass[k] - geometry.cmass[k-1]);
      }

      Pg = (atmos.nHtot[k-1]*atmos.totalAbund + atmos.ne[k-1]) * KBOLTZMANN * atmos.T[k-1];
      Ntot = (Pg + dP) / KBOLTZMANN / atmos.T[k];
      dN = (Ntot - atmos.ne[k]) / atmos.totalAbund - atmos.nHtot[k];
      atmos.nHtot[k] = atmos.nHtot[k]+dN;
      
      if (fabs(dN/atmos.nHtot[k])<break_me)
        break;
    }
    printf("final HSE %d --> %e\n================\n", k, atmos.nHtot[k]);
  }

  for(k=0; k<atmos.Nspace; k++) printf("chi_c[%d] = %e\n", k, chi_c[k]);

  // write_atmos(chi_c, rho);

  return 0;
}

void Solve_ne_(double *ne, bool_t fromscratch, int ki)
{
  const char routineName[] = "Solvene";
  register int k, n, j;

  int     Nmaxstage, niter;
  double *fjk, *dfjk, error, ne_old, akj, sum, PhiH, C1, Uk,
    dne, dnemax, *np, PhiHmin;

  getCPU(3, TIME_START, NULL);

  C1 = (HPLANCK/(2.0*PI*M_ELECTRON)) * (HPLANCK/KBOLTZMANN);

  /* --- Figure out the largest array size needed so that we do not
         have to allocate and free memory all the time -- ----------- */

  Nmaxstage = 0;
  for (n = 0;  n < atmos.Nelem;  n++)
    Nmaxstage = MAX(Nmaxstage, atmos.elements[n].Nstage);
  fjk  = (double *) malloc(Nmaxstage * sizeof(double));
  dfjk = (double *) malloc(Nmaxstage * sizeof(double));

  np = atmos.H->n[atmos.H->Nlevel-1];
  // for (k = 0;  k < atmos.Nspace;  k++) {
  k = ki;
    if (fromscratch) {

      /* --- Get the initial solution from ionization of H only -- -- */

      if (atmos.H_LTE) {
        Uk = getKuruczpf_(&atmos.elements[0], 0, k);
        PhiH = 0.5 * pow(C1/atmos.T[k], 1.5) *
           exp(Uk + atmos.elements[0].ionpot[0]/(KBOLTZMANN*atmos.T[k]));
        ne_old = (sqrt(1.0 + 4.0*atmos.nHtot[k]*PhiH) - 1.0) / (2.0*PhiH);
      } else
         ne_old = np[k];
         /* --- Copy into ne as well to calculate first fij and dfij - - */
         ne[k] = ne_old;
    } else {
      /* --- Use original electron density as starting guess -- ----- */
      ne_old = ne[k];
    }

    niter = 0;
    while (niter < N_MAX_ELECTRON_ITERATIONS) {
      error = ne_old / atmos.nHtot[k];
      sum   = 0.0;
      for (n = 0;  n < atmos.Nelem;  n++) {
        getfjk(&atmos.elements[n], ne_old, k, fjk, dfjk);
        /* --- Contribution from Hminus --             -------------- */
        if (n == 0) {
          PhiHmin = 0.25*pow(C1/atmos.T[k], 1.5) *
      exp(E_ION_HMIN / (KBOLTZMANN * atmos.T[k]));
    error += ne_old * fjk[0] * PhiHmin;
          sum   -= (fjk[0] + ne_old * dfjk[0]) * PhiHmin;
  }
  for (j = 1;  j < atmos.elements[n].Nstage;  j++) {
    akj = atmos.elements[n].abund * j;
    error -= akj * fjk[j];
    sum   += akj * dfjk[j];
  }
      }

      ne[k] = ne_old -
  atmos.nHtot[k] * error / (1.0 - atmos.nHtot[k] * sum);
      dne = fabs((ne[k] - ne_old)/ne_old);
      ne_old = ne[k];
    
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
  // }

  free(fjk);  free(dfjk);

  getCPU(3, TIME_POLL, "Electron density");
}

void backgrOpac_(double *chi_c, double *scatt_c, int ki, double Hmin_fudge, double metal_fudge, double scatt_fudge)
{
  const char routineName[] = "backgrOpac";
  register int n, k;

  double  *chi, *eta, *scatt, *thomson;
  Atom    *He;

  /* --- Only evaluate for static case --              -------------- */

  atmos.moving = FALSE;

  Solve_ne_(atmos.ne, FALSE, ki);
  ChemicalEquilibrium(N_MAX_CHEM_ITER, CHEM_ITER_LIMIT, -1);
  SetLTEQuantities();
  
  chi_c[ki] = 0;
  scatt_c[ki] = 0;

  /* --- Temporary storage for this routine --         -------------- */

  chi   = (double *) malloc(atmos.Nspace * sizeof(double));
  eta   = (double *) malloc(atmos.Nspace * sizeof(double));
  scatt = (double *) malloc(atmos.Nspace * sizeof(double));

  /* --- Check whether He is present among the metals --  ----------- */

  He = (atmos.elements[1].model) ? atmos.elements[1].model : NULL;

  /* --- Thomson scattering by free electrons is wavelength independent
         in non-relativistic limit --                  -------------- */

  thomson = (double *) malloc(atmos.Nspace * sizeof(double));
  Thomson(thomson);
  chi_c[ki] += thomson[ki];

  /* --- Go through wavelengths one-by-one --          -------------- */

  double lambda = atmos.lambda_ref;

	/* --- Negative hydrogen ion, bound-free and free-free -- ------- */

	if (Hminus_bf(lambda, chi, eta)) {
		// for(k=0; k<atmos.Nspace; k++) chi_c[k] += chi[k];
		chi_c[ki] += chi[ki];
	}
	if (Hminus_ff(lambda, chi)) {
		// for(k=0; k<atmos.Nspace; k++) chi_c[k] += chi[k];
    chi_c[ki] += chi[ki];
	}
	
    // for (k = 0;  k < atmos.Nspace;  k++) chi_c[k] *= Hmin_fudge;
  chi_c[ki] *= Hmin_fudge;

	/* --- Bound-free opacities from OH and CH molecules -- --------- */

	// if (OH_bf_opac(lambda, chi, eta)) {
	// 	for(k=0; k<atmos.Nspace; k++) chi_c[k] += chi[k];
	// }
	// if (CH_bf_opac(lambda, chi, eta)) {
	// 	for(k=0; k<atmos.Nspace; k++) chi_c[k] += chi[k];
	// }

	/* --- Neutral Hydrogen Bound-Free and Free-Free --  ------------ */

	if (Hydrogen_bf(lambda, chi, eta)) {
		// for(k=0; k<atmos.Nspace; k++) chi_c[k] += chi[k];
    chi_c[ki] += chi[ki];
	}
	Hydrogen_ff(lambda, chi);
	// for(k=0; k<atmos.Nspace; k++) chi_c[k] += chi[k];
  chi_c[ki] += chi[ki];

	/* --- Rayleigh scattering by neutral hydrogen --  -------------- */

	// if (Rayleigh(lambda, atmos.H, scatt)) {
	// 	for(k=0; k<atmos.Nspace; k++) scatt_c[k] += scatt[k];
	// }
	/* --- Rayleigh scattering by neutral helium --    -------------- */

	// if (He && Rayleigh(lambda, He, scatt)) {
	// 	for(k=0; k<atmos.Nspace; k++) scatt_c[k] += scatt[k];
	// }
	/* --- Absorption by H + H^+ (referred to as H2plus free-free) -- */

	if (H2plus_ff(lambda, chi))
    chi_c[ki] += chi[k];
	// 	for(k=0; k<atmos.Nspace; k++) chi_c[k] += chi[k];

	/* --- Rayleigh scattering and free-free absorption by
	       molecular hydrogen --                       -------------- */

	// if (Rayleigh_H2(lambda, scatt)) {
	// 	for(k=0; k<atmos.Nspace; k++) scatt_c[k] += scatt[k];
	// }
	if (H2minus_ff(lambda, chi)) {
		// for(k=0; k<atmos.Nspace; k++) chi_c[k] += chi[k];
    chi_c[ki] += chi[ki];
	}
	
	/* --- Bound-Free opacities due to ``metals'' --   -------------- */

	// for (n = 1;  n < atmos.Natom;  n++) {
	//   Metal_bf(lambda, 1, &atmos.atoms[n], chi, eta);
	//   for(k=0; k<atmos.Nspace; k++) chi_c[k] += chi[k]*metal_fudge;
	// }

	// for (k = 0;  k < atmos.Nspace;  k++) {
 //      scatt_c[k] *= scatt_fudge;
 //      chi_c[k] += scatt_c[k];
 //    }
    scatt_c[ki] *= scatt_fudge;
    chi_c[ki] += scatt_c[ki];

  /* --- Free the temporary space allocated in the ff routines -- --- */

  // Hminus_ff(0.0, NULL);
  // H2minus_ff(0.0, NULL);
  // H2plus_ff(0.0, NULL);

  free(chi);  free(eta);  free(scatt);//  free(thomson);
}

void write_atmos(double *chi_c, double *rho)
{
	int k;
	char atmos_out[] = "atmosHSE.out";
	double *height, *cmass, *tau_ref;
	FILE *fp_out;
	Atom *H;

	height  = geometry.height;
  	cmass   = geometry.cmass;
  	tau_ref = geometry.tau_ref;

	if (geometry.scale==COLUMN_MASS){
  		// printf("Convert scale\n");
		height[0] = 0.0;
	    cmass[0]  = (tau_ref[0] / chi_c[0]) * rho[0];
	    // printf("%e | %e\n", chi_c[1], rho[1]);
	    for (k = 1;  k < atmos.Nspace;  k++) {
	          height[k] = height[k-1] - 2.0 * (tau_ref[k] - tau_ref[k-1]) / 
	    	(chi_c[k-1] + chi_c[k]);
	          cmass[k]  = cmass[k-1]  + 0.5*(rho[k-1] + rho[k]) *
	    	(height[k-1] - height[k]);
	    }
	}


	fp_out = fopen(atmos_out, "w");

	if (fp_out!=NULL){
		fprintf(fp_out, "* Model file\n");
		fprintf(fp_out, "*\n");
		fprintf(fp_out, "  %s\n", atmos_out);
		fprintf(fp_out, "  Tau scale \n");
		fprintf(fp_out, "*\n* log(g) [cm s^-2]\n");
		fprintf(fp_out, "  %f\n", log10(atmos.gravity/1e2));
		fprintf(fp_out, "*\n* Ndep\n");
		fprintf(fp_out, "  %2d\n", atmos.Nspace);
		fprintf(fp_out, "*\n* log tau    Temp[K]    n_e[cm-3]    v_z[km/s]   v_turb[km/s]\n");
		for(k=0; k<atmos.Nspace; k++){
			fprintf(fp_out, "    %+5.4f    %5.1f   %5.4e   %5.4e   %5.4e\n", log10(geometry.tau_ref[k]), atmos.T[k], atmos.ne[k]/1e6, geometry.vel[k]/1e3, atmos.vturb[k]/1e3);
		}
		fprintf(fp_out, "*\n* Hydrogen populations [cm-3]\n");
		fprintf(fp_out, "*     nh(1)        nh(2)        nh(3)        nh(4)        nh(5)        nh(6)\n");
		for(k=0; k<atmos.Nspace; k++){
			H = &atmos.atoms[0];
			// fprintf(fp_out, "    %5.4e   %5.4e   %5.4e   %5.4e   %5.4e   %5.4e\n", atmos.H->n[0][k]/1e6, atmos.H->n[1][k]/1e6, atmos.H->n[2][k]/1e6, atmos.H->n[3][k]/1e6, atmos.H->n[4][k]/1e6, atmos.H->n[5][k]/1e6);
			fprintf(fp_out, "    %5.4e   %5.4e   %5.4e   %5.4e   %5.4e   %5.4e\n", H->nstar[0][k]/1e6, H->nstar[1][k]/1e6, H->nstar[2][k]/1e6, H->nstar[3][k]/1e6, H->nstar[4][k]/1e6, H->nstar[5][k]/1e6);
		}
	} else{
		printf("Dusan, we have a problem\n");
	}

	fclose(fp_out);
}

double getKuruczpf_(Element *element, int stage, int k)
{
  bool_t hunt = TRUE;
  double Uk;

  Linear(atmos.Npf, atmos.Tpf, element->pf[stage], 
	 1, &atmos.T[k], &Uk, hunt);
  
  return Uk;
}