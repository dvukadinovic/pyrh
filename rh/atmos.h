/* ------- file: -------------------------- atmos.h -----------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Fri Jul  8 15:36:04 2011 --

       --------------------------                      ----------RH-- */

#ifndef __ATMOS_H__
#define __ATMOS_H__

/* --- Define structure to hold geometry-independent
       atmospheric quantities. --                      -------------- */


#define  ATMOS_ID_WIDTH  80

/* --- Maximum values of number of angles per octant in azimuth and
       inclination in case of Gauss-Legendre quadrature in inclination.
       --                                              -------------- */

#define  NMAXINCLINATION  9
#define  NMAXAZIMUTH      5

#define  MOLECULAR_CONCENTRATION_FILE  "molecules.out"

/* --- Angle set identifications. Most are Carlsson type -- --------- */

enum angleset  {SET_VERTICAL, SET_GL, SET_A2, SET_A4, SET_A6, SET_A8,
		SET_B4,	SET_B6, SET_B8, NO_SET};

// enum PERIODIC_SYSTEM_OF_ELEMENTS {"H",  "HE",
//                      "LI", "BE", "B",  "C",  "N",  "O",  "F",  "NE", \
//                      "NA", "MG", "AL", "SI", "P",  "S",  "CL", "AR", \
//                      "K",  "CA", "SC", "TI", "V",  "CR", "MN", "FE", "CO", "NI", "CU", "ZN", \
//                      "GA", "GE", "AS", "SE", "BR", "KR", \
//                      "RB", "SR", "Y",  "ZR", "NB", "MO", "TC", "RU", "RH", "PD", "AG", "CD", \
//                      "IN", "SN", "SB", "TE", "I",  "XE", \
//                      "CS", "BA", \
//                      "LA", "CE", "PR", "ND", "PM", "SM", "EU", "GD", "TB", "DY", "HO", "ER", "TM", "YB", "LU", \
//                      "HF", "TA", "W",  "RE", "OS", "IR", "PT", "AU", "HG", "TL", "PB", "BI", "PO", "AT", "RN", \
//                      "FR", "RA", \
//                      "AC", "TH", "PA", "U", "NP", "PU", "AM", "CM", "BK", "CF", "ES", "FM"};

typedef struct {
  bool_t hasline;
  bool_t ispolarized;
} flags;

typedef struct {
  enum  angleset set;
  int   Ninclination, Nazimuth;
} AngleSet;

typedef struct {
  char    ID[ATMOS_ID_WIDTH];
  bool_t  moving, H_LTE, Stokes, hydrostatic;
  int     Ndim, *N, Nspace, Nrays, Nelem, Natom, Nmolecule,
         *backgrrecno, Npf, NHydr, fd_background, NPRDactive,
          Nactiveatom, Nactivemol;
  int     active_layer; // DV -- layer for which we compute HSE (opacities and equilibrium)
  int     fudge_num; // DV -- number of OF points
  double *fudge_lam, **fudge; // DV -- pointer for OF wavelengths and values
  double *loggf_values, *lam_values; // DV -- perturbation value for log(gf)
  int    *loggf_ids, Nloggf, *lam_ids, Nlam; // DV -- line id for which we perturb log(gf)
  int     Nrlk;
  double *T, *ne, *vturb, totalAbund, avgMolWght, wght_per_H, gravity,
          vmicro_char, vmacro_tresh, lambda_ref, *wmu, *Tpf,
         *nHtot, **nH, *nHmin, *B, *gamma_B, *chi_B, B_char,
        **cos_gamma, **cos_2chi, **sin_2chi;
  double ***atomic_rfs;
  AngleSet  angleSet;
  Element  *elements;
  Atom     *H, *atoms, **activeatoms;
  Molecule *H2, *OH, *CH, *molecules, **activemols;
  RLK_Line *rlk_lines;
  FILE   *fp_atmos;
  flags  *backgrflags;
} Atmosphere;

/* --- Associated function prototypes --               -------------- */

void  freeAtmos(Atmosphere *atmos);
void  readAbundance(Atmosphere *atmos, int Nabun, int *atomic_id, double *atomic_abundance);
void  writeAtmos(Atmosphere *atmos);
void  Solve_ne(double *ne, bool_t fromscratch);
void  initAngleSet(AngleSet *angleSet);


/* --- Background opacities due to various lines --    -------------- */

flags passive_bb(double lambda, int nspect, int mu, bool_t to_obs,
		 double *chi, double *eta, double *chip);

flags rlk_opacity(double lambda, int nspect, int mu, bool_t to_obs,
                  double *chi, double *eta, double *scatt, double *chip);

flags MolecularOpacity(double lambda, int nspect, int mu, bool_t to_obs,
		       double *chi, double *eta, double *chip);


/* --- Moleculear concentrations --                    -------------- */

void writeMolecules(char *fileName);
void readMolecules(char *fileName);


#endif /* !__ATMOS_H__ */

/* ---------------------------------------- atmos.h ----------------- */
