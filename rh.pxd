# cdef extern from "/usr/include/pthread.h":
# ctypedef struct pthread_attr_t_ *pthread_attr_t
# cdef pthread_attr_t dummy

cdef extern from "rh/rh.h":
	cdef double **matrix_double(int Nrow, int Ncol)
	# cdef enum StokesMode: NO_STOKES, FIELD_FREE, POLARIZATION_FREE, FULL_STOKES
	# cdef enum solution: UNKNOWN=-1, LTE_POPULATIONS, ZERO_RADIATION, OLD_POPULATIONS, NEW_J, OLD_J

# cdef extern from "rh/inputs.h":
# 	cdef enum ne_solution: NONE, ONCE, ITERATION
# 	cdef enum order_3D: LINEAR_3D, BICUBIC_3D
# 	cdef enum S_interpol: S_LINEAR, S_PARABOLIC, S_BEZIER3
# 	cdef enum S_interpol_stokes: DELO_PARABOLIC, DELO_BEZIER3

cdef extern from "rh/atom.h":
	cdef enum vdWaals: UNSOLD, RIDDER_RENSBERGEN, BARKLEM, KURUCZ

	ctypedef struct ZeemanMultiplet:
		int Ncomponent
		int *q
		double *shift
		double *strength

	ctypedef struct RLK_Line:
		pass
		# bint polarizable
		# vdWaals vdwaals;
		# int pt_index
		# int stage
		# int isotope
		# int Li
		# int Lj
		# double lambda0, gi, gj, Ei, Ej, Bji, Aji, Bij, Si, Sj
		# double Grad, GStark, GvdWaals, hyperfine_frac
		# double isotope_frac, gL_i, gL_j, hfs_i, hfs_j, iso_dl
		# double cross, alpha
		# ZeemanMultiplet *zm;
	
cdef extern from "rh/rhf1d/pyrh_compute1dray.h":
	ctypedef struct myRLK_Line:
		int Nrlk
		RLK_Line *rlk_lines

	ctypedef struct mySpectrum:
		int nlw
		int Nrays
		int stokes
		double *lam
		double *sI
		double *sQ
		double *sU
		double *sV
		double **J
	
	cdef mySpectrum rhf1d(char *cwd, double mu, int Ndep,
			double *rh_scale, double *rh_temp, double *rh_ne, double *rh_vz, double *rh_vmic,
			double *rh_mag, double *rh_gamma, double *rh_chi,
			double *rh_nH, int atm_scale,
			int Nwave, double *lam,
			int do_fudge, int fudge_num, double *fudge_lam, double *fudge,
			int Nloggf, int *loggf_ids, double* loggf_values,
			int Nlam, int *lam_ids, double *lam_values,
			int NKurucz_lists, char *Kurucz_lists)
	# myRLK_Line *pyrh_rlk_lines,

	cdef myRLK_Line get_RLK_lines(int argc, char *argv[])
	# cdef void dummy(myRLK_Line *rlk_lines)

cdef extern from "rh/rhf1d/pyrh_hse.h":
	ctypedef struct myPops:
		double *ne
		double *nHtot
		double **nH
		double *rho
		double *pg

	cdef myPops hse(char* keyword_input, int Ndep, double pg_top,
					double *rh_scale, double *rh_temp,
					int atm_scale,
					int do_fudge, int fudge_num, double *fudge_lam, double *fudge)

	# cdef void dummy()

	cdef void get_tau(char *cwd, double mu, int pyrh_Ndep, double *tau_ref,
		              double *pyrh_scale, double *pyrh_temp, double *pyrh_ne, double *pyrh_vz, double *pyrh_vmic,
		              double *pyrh_mag, double *pyrh_gamma, double *pyrh_chi,
		              double *pyrh_nH, int pyrh_atm_scale, 
		              double lam_ref)