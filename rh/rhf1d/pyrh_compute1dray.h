#include "../inputs.h"
#include "../atom.h"

typedef struct{
    Atom *atoms;
    Molecule *molecules;
    Element *elements;

    int Nelem; // number of *elements;
    int Npf; // number of partition function elements

    double totalAbund, wght_per_H, avgMolWght;

    double *Tpf; // partition function


} AtMol;

typedef struct{
	int nlw, Nrays, stokes;
	double *lam, *sI, *sQ, *sU, *sV;
	double **J, **J20;
} mySpectrum;

typedef struct{
	int Nrlk;
	RLK_Line *rlk_lines;
} myRLK_Line;

InputData get_InputData(char *cwd);
AtMol read_AtomsMolecules(InputData pyrh_input, char *cwd);
void check_ID(InputData ID);
 
mySpectrum rhf1d(char *cwd, double mu, int Ndep,
              double *pyrh_scale, double *pyrh_temp, double *pyrh_ne, double *pyrh_vz, double *pyrh_vmic,
              double *pyrh_mag, double *pyrh_gamma, double *pyrh_chi,
              double *pyrh_nH, int atm_scale, 
              int Nwave, double *lam,
              int do_fudge, int fudge_num, double *fudge_lam, double *fudge,
              int Nloggf, int *loggf_ids, double *loggf_values,
              int Nlam, int *lam_ids, double *lam_values,
              int NKurucz_lists, char *Kurucz_lists);
              // myRLK_Line *pyrh_rlk_lines,

myRLK_Line get_RLK_lines(int argc, char *argv[]);