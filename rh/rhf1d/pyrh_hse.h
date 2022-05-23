typedef struct{
	double **nH, *ne;
} myPops;

myPops hse(int pyrh_Ndep,
           double *pyrh_scale, double *pyrh_temp, double *pyrh_ne, double *pyrh_vz, double *pyrh_vmic,
           double *pyrh_mag, double *pyrh_gamma, double *pyrh_chi,
           double *pyrh_nH, int pyrh_atm_scale, 
           int do_fudge, int fudge_num, double *fudge_lam, double *fudge);

void dummy();