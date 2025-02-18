/* ------- file: -------------------------- kurucz.c ----------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Thu Oct 12 14:52:27 2017 --

       --------------------------                      ----------RH-- */

/* --- Routines to deal with a Kurucz-format line list. These are lines
       to be used as LTE background.

       The input file should contain a list of Kurucz linelist files,
       one line each, in the format described on the kurucz web page:
         http://kurucz.harvard.edu/linelists.html

FORMAT(F11.4,F7.3,F6.2,F12.3,F5.2,1X,A10,F12.3,F5.2,1X,A10,
3F6.2,A4,2I2,I3,F6.3,I3,F6.3,2I5,1X,A1,A1,1X,A1,A1,i1,A3.2I5,I6) 

 1 wavelength (nm)  air above 200 nm   F11.4
 2 log gf  F7.3
 3 element code = element number + charge/100.  F6.2
 4 first energy level in cm-1   F12.3
 5 J for first level   F5.1
   blank for legibility   1X
 6 label field for first level   A10
 7 second energy level in cm-1   F12.3
        (negative energies are predicted or extrapolated}  
 8 J for second level   F5.1
   blank for legibility   1X
 9 label field for second level   A10
10 log of radiative damping constant, Gamma Rad  F6.2 or F6.3
11 log of stark damping constant/electron number. Gamma Stark  F6.2 or F6.3
12 log of van der Waals damping constant/neutral hydrogen number, 
       Gamma van der Waals   F6.2 or F6.3
13 reference that can be expanded in subdirectory LINES   A4  
14 non-LTE level index for first level   I2
15 non-LTE level index for second level   I2
16 isotope number   I3
17 hyperfine component log fractional strength  F6.3
18 isotope number  (for diatomics there are two and no hyperfine)   I3
19 log isotopic abundance fraction   F6.3
20 hyperfine shift for first level in mK to be added to E  I5
21 hyperfine shift for second level in mK to be added to E'  I5
   the symbol "F" for legibilty   1X
22 hyperfine F for the first level    I1
23 note on character of hyperfine data for first level: z none, ? guessed  A1
   the symbol "-" for legibility    1X
24 hyperfine F' for the second level  I1
25 note on character of hyperfine data for second level: z none, ? guessed  A1
26 1-digit code, sometimes for line strength classes   I1
27 3-character code such as AUT for autoionizing    A3  
28 lande g for first level times 1000   I5
29 lande g for second level times 1000   I5
30 isotope shift of wavelength in mA
31 [DV] (optional) orbital quantum number of the first level I2
32 [DV] (optional) orbital quantum number of the second level I2


 Note: The periodic table index of Kurucz starts at 1 (for Hydrogen).
       In the RH routines counting starts at 0 in the atmos.elements
       structure. This discrepancy is accounted for in the opacity
       calculation (rlk_opacity.c).
 Note: [DV] we added two more values to Kurucz line that contain information
       on the orbital quantum numbers of levels that is used to compute the
       broadening of a spectral line by the ABO theory. The original implementation
       does not read these numbers, it only reads the total angular momentum
       (term state). Do not give these numbers if you do not know the levels
       configuration.

       --                                              -------------- */

#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "background.h"
#include "spectrum.h"
#include "constant.h"
#include "inputs.h"
#include "error.h"

#define COMMENT_CHAR             "#"
// #define RLK_RECORD_LENGTH        160
#define RLK_RECORD_LENGTH        171
#define Q_WING                   20.0
#define MILLI                    1.0E-03
#define ANGSTROM_TO_NM           0.1
#define MAX_GAUSS_DOPPLER        7.0
#define USE_TABULATED_WAVELENGTH 1
#define LN10                     log(10)


/* --- Function prototypes --                          -------------- */

double           RLKProfile(RLK_Line *rlk, int k, int mu, bool_t to_obs,
			    double lambda,
			    double *phi_Q, double *phi_U, double *phi_V,
			    double *psi_Q, double *psi_U, double *psi_V);
ZeemanMultiplet* RLKZeeman(RLK_Line *rlk);
void             initRLK(RLK_Line *rlk);
bool_t           RLKdeterminate(char *labeli, char *labelj, RLK_Line *rlk);
void             getUnsoldcross(RLK_Line *rlk);
void             free_BS(Barklemstruct *bs);
void             freePartitionFunction();


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern InputData input;
extern Spectrum spectrum;
extern char messageStr[];


/* ------- begin -------------------------- readKuruczLines.c ------- */

void readKuruczLines(char *inputFile)
{
  const char routineName[] = "readKuruczLines";
  const double  C = 2.0*PI * (Q_ELECTRON/EPSILON_0) * 
                             (Q_ELECTRON/M_ELECTRON) / CLIGHT;

  char   inputLine[RLK_RECORD_LENGTH+1], listName[MAX_LINE_SIZE],
    filename[MAX_LINE_SIZE], Gvalues[18+1], elem_code[7],
         labeli[RLK_LABEL_LENGTH+1], labelj[RLK_LABEL_LENGTH+1],
        *commentChar = COMMENT_CHAR;
  bool_t swap_levels, determined, useBarklem;
  int    Nline, Nread, Nrequired, checkPoint, hfs_i, hfs_j, gL_i, gL_j,
         iso_dl, line_index;
  double lambda0, Ji, Jj, Grad, GStark, GvdWaals, pti,
         Ei, Ej, gf, lambda_air;
  RLK_Line *rlk;
  Barklemstruct bs_SP, bs_PD, bs_DF;
  FILE  *fp_Kurucz, *fp_linelist;

  if (!strcmp(inputFile, "none")) return;

  /* --- Read in the data files for Barklem collisional broadening -- */

  readBarklemTable(SP, &bs_SP);
  readBarklemTable(PD, &bs_PD);
  readBarklemTable(DF, &bs_DF);

  labeli[RLK_LABEL_LENGTH] = '\0';
  labelj[RLK_LABEL_LENGTH] = '\0';

  if ((fp_Kurucz = fopen(inputFile, "r")) == NULL) {
    sprintf(messageStr, "Unable to open input file %s", inputFile);
    Error(ERROR_LEVEL_1, routineName, messageStr);
    return;
  }

  /* --- Go through each of the linelist files listed in input file - */  

  while (getLine(fp_Kurucz, commentChar, listName, FALSE) != EOF) {
    Nread = sscanf(listName, "%s", filename);
    // update the path to filename?
    if ((fp_linelist = fopen(filename, "r")) == NULL) {
      sprintf(messageStr, "Unable to open input file %s", filename);
      Error(ERROR_LEVEL_1, routineName, messageStr);
    }

    /* --- Count the number of lines in this file --   -------------- */

    Nline = 0;
    while (fgets(inputLine, RLK_RECORD_LENGTH+1, fp_linelist) != NULL)
      if (*inputLine != *commentChar) Nline++;
    rewind(fp_linelist);

    if (atmos.Nrlk == 0) atmos.rlk_lines = NULL;
    atmos.rlk_lines = (RLK_Line *)
      realloc(atmos.rlk_lines, (Nline + atmos.Nrlk) * sizeof(RLK_Line));

    /* --- Read lines from file --                     -------------- */

    rlk = atmos.rlk_lines + atmos.Nrlk;
    line_index = 0;
    while (fgets(inputLine, RLK_RECORD_LENGTH+1, fp_linelist) != NULL) {
      if (*inputLine != *commentChar) {

        initRLK(rlk);

        Nread = sscanf(inputLine, "%lf %lf %s %lf",
      		       &lambda_air, &gf, (char *) &elem_code, &Ei);

        /* --- Ionization stage and periodic table index -- --------- */

        sscanf(elem_code, "%d.%d", &rlk->pt_index, &rlk->stage);

        Nread += sscanf(inputLine+53, "%lf", &Ej);

        Ei = fabs(Ei) * (HPLANCK * CLIGHT) / CM_TO_M;
      	Ej = fabs(Ej) * (HPLANCK * CLIGHT) / CM_TO_M;

        /* --- Beware: the Kurucz linelist has upper and lower levels
	             of a transition in random order. Therefore, we have to
               check for the lowest energy of the two and use that as
               lower level --                          -------------- */

      	if (Ej < Ei) {
      	  swap_levels = TRUE; 
      	  rlk->Ei = Ej;
      	  rlk->Ej = Ei;
      	  strncpy(labeli, inputLine+69, RLK_LABEL_LENGTH);
      	  strncpy(labelj, inputLine+41, RLK_LABEL_LENGTH);
      	} else {
      	  swap_levels = FALSE;
      	  rlk->Ei = Ei;
          rlk->Ej = Ej;
      	  strncpy(labeli, inputLine+41, RLK_LABEL_LENGTH);
      	  strncpy(labelj, inputLine+69, RLK_LABEL_LENGTH);
      	}

      	Nread += sscanf(inputLine+35, "%lf", &Ji);
      	Nread += sscanf(inputLine+63, "%lf", &Jj);
      	if (swap_levels) SWAPDOUBLE(Ji, Jj);
      	rlk->gi = 2*Ji + 1;
      	rlk->gj = 2*Jj + 1;

        if (atmos.Nlam>0){
          for (int idl=0; idl<atmos.Nlam; idl++){
            if (atmos.lam_ids[idl]==line_index){
              lambda_air += atmos.lam_values[idl];
              if (input.get_atomic_rfs){
                rlk->get_dlam_rf = TRUE;
                rlk->dlam_rf_ind = idl;
              }
            }
          }
        }  

      	if (USE_TABULATED_WAVELENGTH) {
      	  /* --- In this case use tabulated wavelength and adjust 
      	         upper-level energy --                 -------------- */

      	  air_to_vacuum(1, &lambda_air, &lambda0);
      	  lambda0 *= NM_TO_M;
      	  rlk->Ej = rlk->Ei +  (HPLANCK * CLIGHT) / lambda0;
      	} else {
      	  /* --- Else use energy levels to calculate lambda0 -- ----- */
      	  
      	  lambda0 = (HPLANCK * CLIGHT) / (rlk->Ej - rlk->Ei);
      	}

        rlk->Aji = C / SQ(lambda0) * POW10(gf) / rlk->gj; 
        if (atmos.Nloggf>=1){
          for (int idl=0; idl<atmos.Nloggf; idl++){
            if (atmos.loggf_ids[idl]==line_index){
              rlk->Aji = C / SQ(lambda0) * POW10(atmos.loggf_values[idl]) / rlk->gj;
              if (input.get_atomic_rfs){
                rlk->get_loggf_rf = TRUE;
                rlk->loggf_rf_ind = idl;
              }
            }
          }
        }
      	rlk->Bji = CUBE(lambda0) / (2.0 * HPLANCK * CLIGHT) * rlk->Aji;
      	rlk->Bij = (rlk->gj / rlk->gi) * rlk->Bji;

        /* --- Store in nm --                          -------------- */

        rlk->lambda0 = lambda0 / NM_TO_M;

        // DV: read orbital quantum numbers of levels 
        int Nread_l, li=-1, lj=-1;
        double alpha, sigma;
        // bool_t got_orbital_numbers = FALSE;
        bool_t got_ABO_coeffs = FALSE;
        Nread_l = sscanf(inputLine+160, "%lf %lf", &alpha, &sigma);
        if (Nread_l!=-1){
          if (((alpha==0.0) || (alpha==1.0) || (alpha==2.0) || (alpha==3.0)) || ((sigma==0.0) || (sigma==1.0) || (sigma==2.0) || (sigma==3.0))){
            // got_orbital_numbers = TRUE;
            li = (int)alpha;
            lj = (int)sigma;
            if (swap_levels){
              SWAPDOUBLE(li, lj);
            }
            rlk->li = li;
            rlk->lj = lj;
          }
          else{
            got_ABO_coeffs = TRUE;
            rlk->alpha = alpha;
            rlk->cross = sigma;
            getABOcross(rlk);
          }
        }

        /* --- Get quantum numbers for angular momentum and spin -- - */

        determined = RLKdeterminate(labeli, labelj, rlk);
        rlk->polarizable = (atmos.Stokes && determined);

        /* --- Line broadening --                      -------------- */

      	strncpy(Gvalues, inputLine+79, 18);
      	Nread += sscanf(Gvalues, "%lf %lf %lf", &Grad, &GStark, &GvdWaals);

      	if (GStark != 0.0) 
      	  rlk->GStark = POW10(GStark) * CUBE(CM_TO_M);
      	else
      	  rlk->GStark = 0.0;

      	if (GvdWaals != 0.0)
      	  rlk->GvdWaals = POW10(GvdWaals) * CUBE(CM_TO_M);
      	else
      	  rlk->GvdWaals = 0.0;

        /* --- If possible use Barklem formalism --    -------------- */

        // useBarklem = FALSE;
        useBarklem = got_ABO_coeffs;

        // [DV] we use orbital quantum number (l) to determine the ABO cross section,
        // not the total angular momentum (L) from the term symbol.
        if (!got_ABO_coeffs){
          if ((rlk->li==0 && rlk->lj==1) || (rlk->li==1 && rlk->lj==0)){
            useBarklem = getBarklemcross(&bs_SP, rlk);
          }
          if ((rlk->li==1 && rlk->lj==2) || (rlk->li==2 && rlk->lj==1)){
            useBarklem = getBarklemcross(&bs_PD, rlk);
          }
          if ((rlk->li==2 && rlk->lj==3) || (rlk->li==3 && rlk->lj==2)){
            useBarklem = getBarklemcross(&bs_DF, rlk);
          }
        }

        if (input.verbose) {
          if (useBarklem)
          {
            printf(" Using ABO broadening for line -- %d @ %f\n", line_index+1, lambda_air);
          } else {
            printf(" No ABO broadening for line    -- %d @ %f\n", line_index+1, lambda_air);
          }
        }

        line_index++;

      	/* --- Else use good old Unsoeld --            -------------- */

        if (!useBarklem) {
      	  getUnsoldcross(rlk);
      	}

      	/* --- Radiative broadening --                 -------------- */

      	if (Grad != 0.0) {
      	  rlk->Grad = POW10(Grad);
      	} else {

      	  /* --- Just take the Einstein Aji value--    -------------- */     

      	  rlk->Grad = rlk->Aji;
      	}

      	/* --- Isotope and hyperfine fractions and slpittings -- ---- */

      	Nread += sscanf(inputLine+106, "%d", &rlk->isotope);
      	Nread += sscanf(inputLine+108, "%lf", &rlk->isotope_frac);
      	rlk->isotope_frac = POW10(rlk->isotope_frac);
      	Nread += sscanf(inputLine+117, "%lf", &rlk->hyperfine_frac);
      	rlk->hyperfine_frac = POW10(rlk->hyperfine_frac);
      	Nread += sscanf(inputLine+123, "%5d%5d", &hfs_i, &hfs_j);
      	rlk->hfs_i = ((double) hfs_i) * MILLI * KBOLTZMANN;
      	rlk->hfs_j = ((double) hfs_j) * MILLI * KBOLTZMANN;

      	/* --- Effective Lande factors --              -------------- */

      	Nread += sscanf(inputLine+143, "%5d%5d", &gL_i, &gL_j);
      	rlk->gL_i = gL_i * MILLI;
      	rlk->gL_j = gL_j * MILLI;
      	if (swap_levels) {
      	  SWAPDOUBLE(rlk->hfs_i, rlk->hfs_j);
      	  SWAPDOUBLE(rlk->gL_i, rlk->gL_j);
      	}

        // [D.V 11.08.2023] Even if SLJ number cannot be read, we can have polarizable line
        //   because user has provided directly the Lande factors for each level.
        //   Works only if LS_LANDE = FALSE (keyword.input file).
        if (!input.LS_Lande && atmos.Stokes){
          if (rlk->gL_i!=-99*MILLI && rlk->gL_j!=-99*MILLI){
            rlk->polarizable = TRUE;
          }
        }

      	/*      Nread += sscanf(inputLine+154, "%d", &iso_dl); */
      	iso_dl = 0;
      	rlk->iso_dl = iso_dl * MILLI * ANGSTROM_TO_NM;

      	checkNread(Nread, Nrequired=17, routineName, checkPoint=1);
	
      	// printf("  Line: %f (vacuum), %f (air)\n"
      	//        " gi, gj: %f, %f\n"
      	//        " Ei, Ej: %e, %e\n"
      	//        " Aji: %e\n"
       //               " Grad, GStark, GvdWaals: %e, %e, %e\n"
       //               " VdWaals: %d\n"
      	//        " hyperfine_frac, isotope_frac: %f, %f\n"
      	//        " cross, alpha: %e, %e\n"
      	//        " Si: %f, Li: %d, Sj: %f, Lj: %d\n"
       //         " gL_i: %f, gL_j: %f\n\n",
      	//        rlk->lambda0, lambda_air,
      	//        rlk->gi, rlk->gj, rlk->Ei, rlk->Ej, rlk->Aji,
      	//        rlk->Grad, rlk->GStark, rlk->GvdWaals,
       //               rlk->vdwaals,
      	//        rlk->hyperfine_frac, rlk->isotope_frac,
      	//        rlk->cross, rlk->alpha, rlk->Si, rlk->Li, rlk->Sj, rlk->Lj,
       //         rlk->gL_i, rlk->gL_j);
      	
      	rlk++;
      } // end of if statement for the text line that does not start with the commentChar
    } // looping through each line in file
    fclose(fp_linelist);

    sprintf(messageStr, "Read %d Kurucz lines from file %s\n",
	    Nline, listName);
    Error(MESSAGE, routineName, messageStr);
    atmos.Nrlk += Nline;
  } // end of looping through list of kurucz linelists

  fclose(fp_Kurucz);

  free_BS(&bs_SP);
  free_BS(&bs_PD);
  free_BS(&bs_DF);
}
/* ------- end ---------------------------- readKuruczLines.c ------- */

/* ------- begin -------------------------- rlk_ascend.c ------------ */
 
int rlk_ascend(const void *v1, const void *v2)
{
  double lambda1 = ((RLK_Line *) v1)->lambda0,
         lambda2 = ((RLK_Line *) v2)->lambda0;

  /* --- Used for sorting transitions by wavelength -- -------------- */

  if (lambda1 < lambda2)
    return -1;
  else if (lambda1 > lambda2)
    return 1;
  else
    return 0;
}
/* ------- end ---------------------------- rlk_ascend.c ------------ */

/* ------- begin -------------------------- rlk_locate.c ------------ */

void rlk_locate(int N, RLK_Line *lines, double lambda, int *low)
{
  int  high, index, increment;

  /* --- Locate position wavelength lambda in Kurucz line list. Assume
         that lines have been sorted in order of ascending wavelength */

  if ((*low <= 0)  ||  (*low > N-1)) {

    /* --- Input guess not useful here, go to bisection --  --------- */

    *low = 0;
    high = N;
  } else {

    /* --- Else hunt up or down to bracket value --    -------------- */ 

    increment = 1;
    if (lambda >= lines[*low].lambda0) {
      high = *low + increment;
      if (*low == N-1) return;

      /* --- Hunt up --                                -------------- */

      while (lambda >= lines[high].lambda0) {
	*low = high;
	increment += increment;
	high = *low + increment;
        if (high >= N) { high = N;  break; }
      }
    } else {
      high = *low;
      if (*low == 0) return;

      /* --- Hunt down --                              -------------- */

      while (lambda <= lines[*low].lambda0) {
	high = *low;
	increment += increment;
	*low = high - increment;
        if (*low <= 0) { *low = 0;  break; }
      }
    }
  }
  /* --- Bisection algorithm --                        -------------- */

  while (high - *low > 1) {
    index = (high + *low) >> 1;
    if (lambda >= lines[index].lambda0)
      *low = index;
    else
      high = index;
  }
}
/* ------- end ---------------------------- rlk_locate.c ------------ */

/* ------- begin -------------------------- rlk_opacity.c ----------- */

flags rlk_opacity(double lambda, int nspect, int mu, bool_t to_obs,
                  double *chi, double *eta, double *scatt, double *chip)
{
  register int k, n, kr;

  bool_t contributes, hunt;
  int    Nwhite, Nblue, Nred, NrecStokes;
  double dlamb_wing, *pf, dlamb_char, hc_la, ni_gi, nj_gj, lambda0, kT,
         Bijhc_4PI, twohnu3_c2, hc, fourPI, hc_4PI,
        *eta_Q, *eta_U, *eta_V, eta_l,
        *chi_Q, *chi_U, *chi_V, chi_l, *chip_Q, *chip_U, *chip_V,
         phi, phi_Q, phi_U, phi_V, psi_Q, psi_U, psi_V,
         epsilon, C, C2_atom, C2_ion, C3, dE, x;
  Atom *metal;
  AtomicLine *line;
  Element *element;
  RLK_Line *rlk;
  flags backgrflags;

  /* --- Calculate the LTE opacity at wavelength lambda due to atomic
         transitions stored in atmos.rlk_lines --      -------------- */

  backgrflags.hasline     = FALSE;
  backgrflags.ispolarized = FALSE;

  /* --- If wavelength outside our list return without calculation -- */

  dlamb_char = lambda * Q_WING * (atmos.vmicro_char / CLIGHT);
  if (lambda < atmos.rlk_lines[0].lambda0 - dlamb_char ||
      lambda > atmos.rlk_lines[atmos.Nrlk-1].lambda0 + dlamb_char) {
   return backgrflags;
  }

  hc     = HPLANCK * CLIGHT;
  fourPI = 4.0 * PI;
  hc_4PI = hc / fourPI;

  if (input.rlkscatter) {
    C       = 2 * PI * (Q_ELECTRON/EPSILON_0) *
                (Q_ELECTRON/M_ELECTRON) / CLIGHT;
    C2_atom = 2.15E-6;
    C2_ion  = 3.96E-6;
  }

  pf = (double *) malloc(atmos.Nspace * sizeof(double));

  /* --- locate wavelength lambda in table of lines -- -------------- */

  Nwhite = 0;
  rlk_locate(atmos.Nrlk, atmos.rlk_lines, lambda, &Nwhite);
  Nblue = Nwhite;
  while (atmos.rlk_lines[Nblue].lambda0 + dlamb_char > lambda &&
	 Nblue > 0)  Nblue--;
  Nred = Nwhite;
  while (atmos.rlk_lines[Nred].lambda0 - dlamb_char < lambda &&
	 Nred < atmos.Nrlk-1)  Nred++;

  /* --- Initialize the contribution for this wavelength and angle -- */

  if (Nred >= Nblue) {
    if (atmos.Stokes) {
      NrecStokes = 4;

      /* --- Use pointers to sub-arrays for Q, U, and V -- ---------- */

      chi_Q = chi + atmos.Nspace;
      chi_U = chi + 2*atmos.Nspace;
      chi_V = chi + 3*atmos.Nspace;

      eta_Q = eta + atmos.Nspace;
      eta_U = eta + 2*atmos.Nspace;
      eta_V = eta + 3*atmos.Nspace;

      if (input.magneto_optical) {
        chip_Q = chip;
        chip_U = chip + atmos.Nspace;
        chip_V = chip + 2*atmos.Nspace;

        for (k = 0;  k < 3*atmos.Nspace;  k++) chip[k] = 0.0;
      }
    } else
      NrecStokes = 1;

    for (k = 0;  k < NrecStokes * atmos.Nspace;  k++) {
      chi[k] = 0.0;
      eta[k] = 0.0;
    }
    if (input.rlkscatter) {
      for (k = 0;  k < atmos.Nspace;  k++) scatt[k] = 0.0;
    }
  }

  /* --- Add opacities from lines at this wavelength -- ------------- */

  for (n = Nblue;  n <= Nred;  n++) {
    rlk = &atmos.rlk_lines[n];
    
    if (fabs(rlk->lambda0 - lambda) <= dlamb_char) {      
      element = &atmos.elements[rlk->pt_index - 1];

      /* --- Check whether partition function is present for this
	     stage, and if abundance is set --         -------------- */

      if ((rlk->stage < element->Nstage - 1) && element->abundance_set) {
      	contributes = TRUE;
      	if ((metal = element->model) != NULL) {

          /* --- If an explicit atomic model is present check that we
	        do not already account for this line in this way - - */

      	  for (kr = 0;  kr < metal->Nline;  kr++) {
      	    line = metal->line + kr;
      	    dlamb_wing = line->lambda0 * line->qwing *
      	      (atmos.vmicro_char / CLIGHT);
      	    if (fabs(lambda - line->lambda0) <= dlamb_wing &&
      		    metal->stage[line->i] == rlk->stage) {
      	      contributes = FALSE;
      	      break;
      	    }
      	  }
      	}
      } else
        contributes = FALSE;

      /* --- Get opacity from line --                  -------------- */

      if (contributes) {
      	hc_la      = (HPLANCK * CLIGHT) / (rlk->lambda0 * NM_TO_M);
      	Bijhc_4PI  = hc_4PI * rlk->Bij * rlk->isotope_frac * rlk->hyperfine_frac * rlk->gi;
      	twohnu3_c2 = rlk->Aji / rlk->Bji;

      	if (input.rlkscatter) {
      	  if (rlk->stage == 0) {
      	    x  = 0.68;
      	    C3 = C / (C2_atom * SQ(rlk->lambda0 * NM_TO_M));
      	  } else {
      	    x  = 0.0;
      	    C3 = C / (C2_ion * SQ(rlk->lambda0 * NM_TO_M));
      	  }

      	  dE = rlk->Ej - rlk->Ei;
      	}
        
        /* --- Set flag that line is present at this wavelength -- -- */

      	backgrflags.hasline = TRUE;
      	if (rlk->polarizable) {
      	  backgrflags.ispolarized = TRUE;
      	  if (rlk->zm == NULL) rlk->zm = RLKZeeman(rlk);
      	}

        if (element->n == NULL) {
      	  element->n = matrix_double(element->Nstage, atmos.Nspace);
      	  LTEpops_elem(element);
      	}
        Linear(atmos.Npf, atmos.Tpf, element->pf[rlk->stage],
         atmos.Nspace, atmos.T, pf, hunt=TRUE);

      	for (k = 0;  k < atmos.Nspace;  k++) {
      	  phi = RLKProfile(rlk, k, mu, to_obs, lambda,
      			   &phi_Q, &phi_U, &phi_V,
      			   &psi_Q, &psi_U, &psi_V);

      	  if (phi){
      	    kT    = 1.0 / (KBOLTZMANN * atmos.T[k]);
      	    ni_gi = element->n[rlk->stage][k] * exp(-rlk->Ei*kT - pf[k]);
            nj_gj = ni_gi * exp(-hc_la * kT);

      	    chi_l = Bijhc_4PI * (ni_gi - nj_gj);
      	    eta_l = Bijhc_4PI * twohnu3_c2 * nj_gj;

      	    if (input.rlkscatter) {
      	      epsilon = 1.0 / (1.0 + C3 * pow(atmos.T[k], 1.5) /
      			       (atmos.ne[k] * pow(KBOLTZMANN * atmos.T[k] / dE, 1 + x)));

              scatt[k] += (1.0 - epsilon) * chi_l * phi;
      	      chi_l    *= epsilon; 
              eta_l    *= epsilon;

              // if (rlk->get_loggf_rf) dscatt[k][rlk->loggf_rf_ind] = scatt[k] * LN10; // this was done from head, for log(gf) should be correct
      	    }

      	    chi[k] += chi_l * phi;
      	    eta[k] += eta_l * phi;

            if (rlk->get_loggf_rf){
              spectrum.dchi_c_lam[nspect][k][rlk->loggf_rf_ind] = chi_l * phi * LN10;
              spectrum.deta_c_lam[nspect][k][rlk->loggf_rf_ind] = eta_l * phi * LN10;
            }

      	    if (rlk->zm != NULL && rlk->Grad) {
      	      chi_Q[k] += chi_l * phi_Q;
      	      chi_U[k] += chi_l * phi_U;
      	      chi_V[k] += chi_l * phi_V;

      	      eta_Q[k] += eta_l * phi_Q;
      	      eta_U[k] += eta_l * phi_U;
      	      eta_V[k] += eta_l * phi_V;

      	      if (input.magneto_optical) {
            		chip_Q[k] += chi_l * psi_Q;
            		chip_U[k] += chi_l * psi_U;
            		chip_V[k] += chi_l * psi_V;
      	      }
      	    }
      	  } // end if(phi)
      	} // end for loop over k
      } // end if(contributes)
    } // end if (lam <= dlamchar)
  } // end for loop over n

  free(pf);
  return backgrflags;
}
/* ------- end ---------------------------- rlk_opacity.c ----------- */

/* ------- begin -------------------------- RLKProfile.c ------------ */

double RLKProfile(RLK_Line *rlk, int k, int mu, bool_t to_obs,
                  double lambda,
		  double *phi_Q, double *phi_U, double *phi_V,
		  double *psi_Q, double *psi_U, double *psi_V)
{
  register int nz;

  double v, phi_sm, phi_sp, phi_pi, psi_sm, psi_sp, psi_pi, adamp,
         vB, H, F, sv, phi_sigma, phi_delta, sign, sin2_gamma, phi,
         psi_sigma, psi_delta, vbroad, vtherm, GvdW, *np;
  Element *element;

  /* --- Returns the normalized profile for a Kurucz line
         and calculates the Stokes profile components if necessary -- */

  element = &atmos.elements[rlk->pt_index - 1];
  vtherm  = 2.0*KBOLTZMANN/(AMU * element->weight);
  vbroad  = sqrt(vtherm*atmos.T[k] + SQ(atmos.vturb[k]));

  v = (lambda/rlk->lambda0 - 1.0) * CLIGHT/vbroad;
  if (atmos.moving) {
    if (to_obs)
      v += vproject(k, mu) / vbroad;
    else
      v -= vproject(k, mu) / vbroad;
  }

  sv = 1.0 / (SQRTPI * vbroad);

  if (rlk->Grad) {
    switch (rlk->vdwaals) {
    case UNSOLD:
      GvdW = rlk->cross * pow(atmos.T[k], 0.3);
      break;

    case BARKLEM:
      GvdW = rlk->cross * pow(atmos.T[k], (1.0 - rlk->alpha)/2.0);
      break;

    default:
      GvdW = rlk->GvdWaals;
      break;
    }
    np = atmos.H->n[atmos.H->Nlevel-1];
    adamp = (rlk->Grad + rlk->GStark * atmos.ne[k] + 
	     GvdW * (atmos.nHtot[k] - np[k])) * 
      (rlk->lambda0  * NM_TO_M) / (4.0*PI * vbroad);
  } else {
    phi = (fabs(v) <= MAX_GAUSS_DOPPLER) ? exp(-v*v) : 0.0;
    return phi * sv;
  }

  if (rlk->polarizable) {
    sin2_gamma = 1.0 - SQ(atmos.cos_gamma[mu][k]);
    vB   = (LARMOR * rlk->lambda0) * atmos.B[k] / vbroad;
    sign = (to_obs) ? 1.0 : -1.0;

    phi_sm = phi_pi = phi_sp = 0.0;
    psi_sm = psi_pi = psi_sp = 0.0;

    for (nz = 0;  nz < rlk->zm->Ncomponent;  nz++) {
      H = Voigt(adamp, v - rlk->zm->shift[nz]*vB, &F, HUMLICEK);

      switch (rlk->zm->q[nz]) {
      case -1:
      	phi_sm += rlk->zm->strength[nz] * H;
      	psi_sm += rlk->zm->strength[nz] * F;
      	break;
      case  0:
      	phi_pi += rlk->zm->strength[nz] * H;
      	psi_pi += rlk->zm->strength[nz] * F;
      	break;
      case  1:
      	phi_sp += rlk->zm->strength[nz] * H;
      	psi_sp += rlk->zm->strength[nz] * F;
      }
    }
    phi_sigma = phi_sp + phi_sm;
    phi_delta = 0.5*phi_pi - 0.25*phi_sigma;

    phi = (phi_delta*sin2_gamma + 0.5*phi_sigma) * sv;

    *phi_Q = sign * phi_delta * sin2_gamma * atmos.cos_2chi[mu][k] * sv;
    *phi_U = phi_delta * sin2_gamma * atmos.sin_2chi[mu][k] * sv;
    *phi_V = sign * 0.5*(phi_sp - phi_sm) * atmos.cos_gamma[mu][k] * sv;

    if (input.magneto_optical) {
      psi_sigma = psi_sp + psi_sm;
      psi_delta = 0.5*psi_pi - 0.25*psi_sigma;

      *psi_Q = sign * psi_delta * sin2_gamma * atmos.cos_2chi[mu][k] * sv;
      *psi_U = psi_delta * sin2_gamma * atmos.sin_2chi[mu][k] * sv;
      *psi_V = sign * 0.5*(psi_sp - psi_sm) * atmos.cos_gamma[mu][k] * sv;
    }
  } else
   phi = Voigt(adamp, v, NULL, ARMSTRONG) * sv;

  return phi;
}
/* ------- end ---------------------------- RLKProfile.c ------------ */

/* ------- begin -------------------------- RLKZeeman.c ------------- */

ZeemanMultiplet* RLKZeeman(RLK_Line *rlk)
{
  const char routineName[] = "RLKZeeman";

  register int n;

  double Jl, Ju, Mu, Ml, norm[3], gLu, gLl, g_eff;
  ZeemanMultiplet *zm;

  /* --- Return a pointer to a ZeemanMultiplet structure with all the
         components of a Zeeman split line. The strengths in the line
         are normalized to unity for each of the three possible values
         of q = [-1, 0, 1].

         Convention:
 
          -- q = +1 corresponds to a redshifted \sigma profile
	     (zm->shift > 0). This redshifted profile has
             right-handed circular polarization when the
             magnetic field parallel to the line of sight and
             points towards the observer.

          -- q = 0 corresponds to an unpolarized \pi profile

	 --                                            -------------- */

  Jl = (rlk->gi - 1.0) / 2.0;
  Ju = (rlk->gj - 1.0) / 2.0;
  zm = (ZeemanMultiplet *) malloc(sizeof(ZeemanMultiplet));

  /* --- Count the number of components --           -------------- */

  zm->Ncomponent = 0;
  for (Ml = -Jl;  Ml <= Jl;  Ml++) {
    for (Mu = -Ju;  Mu <= Ju;  Mu++)
      if (fabs(Mu - Ml) <= 1.0) zm->Ncomponent++;
  }

  zm->q        = (int *) malloc(zm->Ncomponent * sizeof(int));
  zm->strength = (double *) malloc(zm->Ncomponent * sizeof(double));
  zm->shift    = (double *) malloc(zm->Ncomponent * sizeof(double));

  for (n = 0;  n < 3;  n++) norm[n] = 0.0;
  g_eff = 0.0;

  /* --- Fill the structure and normalize the strengths -- -------- */

  if (input.LS_Lande){
    /* Lande factors in LS coupling */
    gLl = Lande(rlk->Si, rlk->Li, Jl);
    gLu = Lande(rlk->Sj, rlk->Lj, Ju);
  }
  else{
    /* Compute Lande factors for each level in the LS coupling 
       if -99 values are provided in the Kurucz line list.
    */
    if (rlk->gL_i==-99*MILLI){
      gLl = Lande(rlk->Si, rlk->Li, Jl);
    } else {
      gLl = rlk->gL_i;
    }

    if (rlk->gL_j==-99*MILLI){
      gLu = Lande(rlk->Sj, rlk->Lj, Ju);
    } else {
      gLu = rlk->gL_j;
    }
  }

  n = 0;
  for (Ml = -Jl;  Ml <= Jl;  Ml++) {
    for (Mu = -Ju;  Mu <= Ju;  Mu++) {
      if (fabs(Mu - Ml) <= 1.0) {
      	zm->q[n]        = (int) (Ml - Mu);
      	zm->shift[n]    = gLl*Ml - gLu*Mu;
      	zm->strength[n] = ZeemanStrength(Ju, Mu, Jl, Ml);
      	  
      	norm[zm->q[n]+1] += zm->strength[n];
        if (zm->q[n] == 1) g_eff += zm->shift[n] * zm->strength[n];
      	n++;
      }
    }
  }
  for (n = 0;  n < zm->Ncomponent;  n++)
    zm->strength[n] /= norm[zm->q[n]+1];
  g_eff /= norm[2];

  return zm;
}
/* ------- end ---------------------------- RLKZeeman.c ------------- */

/* ------- begin -------------------------- RLKdeterminate.c -------- */

bool_t RLKdeterminate(char *labeli, char *labelj, RLK_Line *rlk)
{
  const char routineName[] = "RLKZeeman";

  char **words, orbit[2];
  bool_t invalid;
  int    count, multiplicity, length, Nread, Ji, Jj;

  /* --- Get spin and orbital quantum numbers from level labels -- -- */

  words  = getWords(labeli, " ", &count);
  if (words[0]) {
    length = strlen(words[count-1]);
    Nread  = sscanf(words[count-1] + length-2, "%d%1s",
		    &multiplicity, orbit);
    free(words);
    if (Nread != 2 || !isupper(orbit[0])) return FALSE;
    
    rlk->Li = getOrbital(orbit[0]);
    rlk->Si = (multiplicity - 1) / 2.0;
    Ji = (rlk->gi - 1.0) / 2.0;
  } else
    return FALSE;

  words  = getWords(labelj, " ", &count);
  if (words[0]) {
    length = strlen(words[count-1]);
    Nread  = sscanf(words[count-1] + length-2, "%d%1s",
		    &multiplicity, orbit);
    free(words);
    if (Nread != 2 || !isupper(orbit[0])) return FALSE;

    rlk->Lj = getOrbital(orbit[0]);
    rlk->Sj = (multiplicity - 1) / 2.0;
    Jj = (rlk->gj - 1.0) / 2.0;
  } else
    return FALSE;

  /* --- For the moment only allow electronic dipole transitions -- --*/ 

  /*  if (fabs(Ji - Jj) > 1.0)
    return FALSE;
    else */
    return TRUE;
}
/* ------- end ---------------------------- RLKdeterminate.c -------- */

/* ------- begin -------------------------- initRLK.c --------------- */

void initRLK(RLK_Line *rlk)
{
  rlk->polarizable = FALSE;
  rlk->zm = NULL; 
  rlk->get_loggf_rf = FALSE;
  rlk->get_dlam_rf = FALSE;
}
/* ------- end ---------------------------- initRLK.c --------------- */

/* ------- begin -------------------------- getUnsoldcross.c -------- */

void getUnsoldcross(RLK_Line *rlk)
{
  const double FOURPIEPS0 = 4.0 * PI * EPSILON_0;

  double   Z, deltaR, vrel35_H, vrel35_He, C625;
  Element *element, *H, *He;

  element = &atmos.elements[rlk->pt_index - 1];
  H = &atmos.elements[0];
  He = &atmos.elements[1];

  if (rlk->stage > element->Nstage - 1) {
    rlk->vdwaals = KURUCZ;
    return;
  }
  
  Z = rlk->stage + 1;
  deltaR = SQ(E_RYDBERG/(element->ionpot[rlk->stage] - rlk->Ej)) -
    SQ(E_RYDBERG/(element->ionpot[rlk->stage] - rlk->Ei));

  if (deltaR <= 0.0) {
    rlk->vdwaals = KURUCZ;
    return;
  }

  vrel35_H  = pow(8.0*KBOLTZMANN/(PI * AMU * element->weight) * 
		  (1.0 + element->weight/H->weight), 0.3);
  vrel35_He = pow(8.0*KBOLTZMANN/(PI * AMU * element->weight) * 
		  (1.0 + element->weight/He->weight), 0.3);
  C625 = pow(2.5 * (SQ(Q_ELECTRON)/FOURPIEPS0) *
	     (ABARH/FOURPIEPS0) *
	     2*PI * SQ(Z*RBOHR)/HPLANCK * deltaR, 0.4);
  rlk->cross = 8.08 *(vrel35_H + He->abund*vrel35_He) * C625;

  rlk->vdwaals = UNSOLD;
}
/* ------- end ---------------------------- getUnsoldcross.c -------- */

/* ------- begin -------------------------- free_BS.c --------------- */

void free_BS(Barklemstruct *bs)
{
  free(bs->neff1);
  free(bs->neff2);

  freeMatrix((void **) bs->cross);
  freeMatrix((void **) bs->alpha);
}
/* ------- end ---------------------------- free_BS.c --------------- */

void freePartitionFunction()
{
  free(atmos.Tpf);  atmos.Tpf = NULL;
  for (int n = 0;  n < atmos.Nelem;  n++) {
    free(atmos.elements[n].mol_index);
    free(atmos.elements[n].ionpot);
    freeMatrix((void **) atmos.elements[n].pf);
    if (atmos.elements[n].n) freeMatrix((void **) atmos.elements[n].n);
  }
}