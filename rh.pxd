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
	
cdef extern from "rh/rhf1d/pyrh_compute1dray.h":
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
	
	cdef mySpectrum rhf1d(int argc, char *argv[], int Ndep,
			double *rh_scale, double *rh_temp, double *rh_ne, double *rh_vz, double *rh_vmic,
			double *rh_mag, double *rh_gamma, double *rh_chi,
			double **rh_nH,
			int atm_scale)

	cdef void get_RLK_lines(int argc, char *argv[])