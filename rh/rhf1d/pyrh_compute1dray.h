#include "../inputs.h"
#include "../atom.h"

typedef struct{
    char ID[10];
    int Nlevel;
    int Nz;
    double **n;
    double **nstar;
} AtomPops;

typedef struct{
	int nlw, Nrays, stokes;
	double *lam, *sI, *sQ, *sU, *sV;
	double **J, **J20;
    double **rfs;
    int Nactive_atoms;
    AtomPops *atom_pops;
} mySpectrum;

typedef struct{
	int Nrlk;
	RLK_Line *rlk_lines;
} myRLK_Line;
 
mySpectrum rhf1d(char *cwd, double mu, int Ndep,
              double *pyrh_scale, double *pyrh_temp, double *pyrh_ne, double *pyrh_vz, double *pyrh_vmic,
              double *pyrh_mag, double *pyrh_gamma, double *pyrh_chi,
              double *pyrh_nH, int atm_scale, 
              int Nwave, double *lam,
              int fudge_num, double *fudge_lam, double *fudge,
              int Nloggf, int *loggf_ids, double *loggf_values,
              int Nlam, int *lam_ids, double *lam_values,
              int Nabun, int *atomic_id, double *atomic_abundance,
              int get_atomic_rfs, int get_populations,
              int NKurucz_lists, char *Kurucz_lists);
              // myRLK_Line *pyrh_rlk_lines,

myRLK_Line get_RLK_lines(char *cwd);
