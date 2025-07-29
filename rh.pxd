# cdef extern from "/usr/include/pthread.h":
# ctypedef struct pthread_attr_t_ *pthread_attr_t
# cdef pthread_attr_t dummy

ctypedef struct Dusan:
	int n
	double dusan
	double *run
	# int **again

cdef extern from "headers/types.h":
	ctypedef int bool_t

cdef extern from "rh/rh.h":
	cdef double **matrix_double(int Nrow, int Ncol)
	cdef int **matrix_int(int Nrow, int Ncol)
	cdef void freeMatrix(void **matrix)
	cdef enum StokesMode: NO_STOKES, FIELD_FREE, POLARIZATION_FREE, FULL_STOKES
	cdef enum solution: UNKNOWN=-1, LTE_POPULATIONS, ZERO_RADIATION, OLD_POPULATIONS, NEW_J, OLD_J
	cdef double** matrix_double(int Nrow, int Ncol)

cdef extern from "rh/inputs.h":
	cdef enum ne_solution: NONE, ONCE, ITERATION
	cdef enum S_interpol: S_LINEAR, S_PARABOLIC, S_BEZIER3
	cdef enum S_interpol_stokes: DELO_PARABOLIC, DELO_BEZIER3

# produced error because geometry.h does not know Atmosphere structure...
# cdef extern from "rh/rhf1d/geometry.h":
# 	cdef enum mass_scale: GEOMETRIC, COLUMN_MASS, TAU500

cdef extern from "rh/atom.h":
	ctypedef struct Atom:
		pass

	ctypedef struct Element:
		int Nstage
		double *ionpot
		double **pf

	ctypedef struct Molecule:
		pass

	cdef enum vdWaals: UNSOLD, RIDDER_RENSBERGEN, BARKLEM, KURUCZ

	ctypedef struct ZeemanMultiplet:
		int Ncomponent
		int *q
		double *shift
		double *strength

	ctypedef struct RLK_Line:
		pass
		#int loggf_rf_ind;
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

cdef extern from "rh/atmos.h":
	ctypedef struct Atmosphere:
		double totalAbund, wght_per_H, avgMolWght
		int Npf, Nelem
		double *Tpf
		Element *elements

cdef extern from "rh/inputs.h":
	ctypedef struct InputData:
		int isum, Ngdelay, Ngorder, Ngperiod, NmaxIter
		int PRD_NmaxIter, PRD_Ngdelay, PRD_Ngorder, PRD_Ngperiod
		int NmaxScatter, Nthreads
		int n_atomic_pars # number of atomic line parameters for which we want to compute RFs
		double iterLimit, PRDiterLimit, metallicity
		double *abundances
		ne_solution solve_ne
		S_interpol_stokes S_interpolation_stokes
		S_interpol S_interpolation
		StokesMode StokesMode
		solution startJ
		bool_t magneto_optical
		bool_t PRD_angle_dep
		bool_t XRD
		bool_t Eddington
		bool_t backgr_pol
		bool_t limit_memory
		bool_t allow_passive_bb
		bool_t NonICE
		bool_t rlkscatter
		bool_t xdr_endian
		bool_t old_background
		bool_t accelerate_mols
		bool_t do_fudge
		bool_t pyrhHSE
		bool_t get_atomic_rfs
		bool_t LS_Lande
		bool_t solve_NLTE
		bool_t verbose
		bool_t get_populations
		bool_t read_atom_model
		char abund_input[160]

cdef extern from "rh/rhf1d/pyrh_compute1dray.h":
	ctypedef struct myRLK_Line:
		int Nrlk
		RLK_Line *rlk_lines

	ctypedef struct AtomPops:
		char ID[10]
		int Nlevel
		int Nz
		double **n
		double **nstar

	ctypedef struct mySpectrum:
		int nlw
		int Nrays
		int stokes
		int Nactive_atoms
		double *lam
		double *sI
		double *sQ
		double *sU
		double *sV
		double **J
		double **rfs
		AtomPops *atom_pops

	ctypedef struct AtMol:
		pass
		# Atom *atoms;
		# Molecule *molecules;
		# Element *elements;
	
	cdef mySpectrum rhf1d(char *cwd, double mu, int Ndep,
			double *rh_scale, double *rh_temp, double *rh_ne, double *rh_vz, double *rh_vmic,
			double *rh_mag, double *rh_gamma, double *rh_chi,
			double *rh_nH, int atm_scale,
			int Nwave, double *lam,
			int fudge_num, double *fudge_lam, double *fudge,
			int Nloggf, int *loggf_ids, double* loggf_values,
			int Nlam, int *lam_ids, double *lam_values,
			int Nabun, int *abundance_id, double *abundance_value,
			int get_atomic_rfs, int get_populations,
			int NKurucz_lists, char *Kurucz_lists)
	# myRLK_Line *pyrh_rlk_lines,

	cdef myRLK_Line get_RLK_lines(char *cwd)
	# cdef void dummy(myRLK_Line *rlk_lines)
	cdef void read_inputs(char *cwd, InputData *pyrh_input_data, Atmosphere *pyrh_atmos)
	cdef void check_inputs(InputData ID, Atmosphere a)

cdef extern from "rh/rhf1d/pyrh_hse.h":
	ctypedef struct myPops:
		double *ne
		double *nHtot
		double **nH
		double *rho
		double *pg

	cdef myPops hse(char* keyword_input, int Ndep,
					double *rh_scale, double *rh_temp,
					double *rh_ne, double *rh_nHtot, double *rh_rho, double *rh_pg, 
					int atm_scale,
					int fudge_num, double *fudge_lam, double *fudge,
					int Nabun, int *abundance_id, double *abundance_value)

	# cdef void dummy()

	cdef void get_scales(char *cwd, int pyrh_Ndep,
			              double *pyrh_scale, double *pyrh_temp, double *pyrh_ne, double *pyrh_vz, double *pyrh_vmic,
			              double *pyrh_nH, int pyrh_atm_scale, 
			              double lam_ref, double *tau, double *height, double *cmass,
						  int Nabun, int *abundance_id, double *abundance_value,)

	cdef void get_ne_from_nH(char *cwd, 
					int pyrh_atm_scale, int pyrh_Ndep, 
                    double *pyrh_scale, double *pyrh_temp, 
                    double *pyrh_nH, double *pyrh_ne)

cdef extern from "rh/rhf1d/pyrh_read_input.h":
	cdef void test_InputData(InputData pyrh_input, Atmosphere pyrh_atmos)
	cdef void set_elements(InputData pyrh_input, Atmosphere *pyrh_atmos)
	cdef void read_atom_model(char *cwd, char *filename)
