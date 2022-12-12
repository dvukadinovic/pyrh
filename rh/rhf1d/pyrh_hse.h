typedef struct{
	double **nH, *ne, *nHtot, *rho, *pg;
} myPops;

myPops hse(char* cwd, int pyrh_Ndep, double pg_top,
           double *pyrh_scale, double *pyrh_temp, 
           int pyrh_atm_scale, 
           int do_fudge, int fudge_num, double *fudge_lam, double *fudge);

void dummy();