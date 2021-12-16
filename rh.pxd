# cdef extern from "/usr/include/pthread.h":
# ctypedef struct pthread_attr_t_ *pthread_attr_t
# cdef pthread_attr_t dummy

cdef extern from "rh/rh.h":
	cdef double **matrix_double(int Nrow, int Ncol)
	cdef enum StokesMode: NO_STOKES, FIELD_FREE, POLARIZATION_FREE, FULL_STOKES
	cdef enum solution: UNKNOWN=-1, LTE_POPULATIONS, ZERO_RADIATION, OLD_POPULATIONS, NEW_J, OLD_J

cdef extern from "rh/inputs.h":
	cdef enum ne_solution: NONE, ONCE, ITERATION
	cdef enum order_3D: LINEAR_3D, BICUBIC_3D
	cdef enum S_interpol: S_LINEAR, S_PARABOLIC, S_BEZIER3
	cdef enum S_interpol_stokes: DELO_PARABOLIC, DELO_BEZIER3
	
	ctypedef struct InputData:
		char atmos_input[160]
		char abund_input[160]
		char wavetable_input[160]
		char atoms_input[160]
		char molecules_input[160]
		char Stokes_input[160]
		char KuruczData[160]
		char pfData[160]
		char fudgeData[160]
		char atmos_output[160]
		char spectrum_output[160]
		char geometry_output[160]
		char opac_output[160]
		char JFile[160]
		char background_File[160]
		char background_ray_File[160]
		char H_atom[160]
		char H2_molecule[160]
		char radrateFile[160]
		char collrateFile[160]
		char dampingFile[160]
		char coolingFile[160]
		char Itop[160]
		bint magneto_optical
		bint PRD_angle_dep
		bint XRD, Eddington
		bint backgr_pol
		bint limit_memory
		bint allow_passive_bb 
		bint NonICE
		bint rlkscatter
		bint xdr_endian
		bint old_background
		bint accelerate_mols,
		solution startJ;
	  
		StokesMode StokesMode;
		S_interpol S_interpolation;
		S_interpol_stokes S_interpolation_stokes;
	  
		order_3D interpolate_3D
		ne_solution solve_ne
		int isum, Ngdelay, Ngorder, Ngperiod, NmaxIter
		int	PRD_NmaxIter, PRD_Ngdelay, PRD_Ngorder, PRD_Ngperiod
		int	NmaxScatter, Nthreads
		double iterLimit, PRDiterLimit, metallicity;

		# pthread_attr_t thread_attr;

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
	
	cdef InputData readMe(int argc, char *argv[])