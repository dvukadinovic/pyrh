# distutils: language = c
# distutils: sources = [rh/rhf1d/pyrh_compute1dray.c]
# distutils: include_dirs = rh/

__version__ = 0.30

cimport rh
import numpy as np
cimport numpy as cnp
import cython
import ctypes
import sys
import time
import xdrlib

from libc.string cimport memcpy
from libc.stdlib cimport malloc, free

cnp.import_array()

default_abundances = np.array([12.00, 10.99, 1.16, 1.15, 2.60, 8.39, 8.00, # H
					    8.66,  4.40, 8.09, 6.33, 7.58, 6.47, 7.55, # O...
						5.45,  7.21, 5.50, 6.56, 5.12, 6.36, 3.10, # P...
						4.99,  4.00, 5.67, 5.39, 7.44, 4.92, 6.25, # Ti...
						4.21,  4.60, 2.88, 3.41, 2.37, 3.35, 2.63, # Cu...
						3.23,  2.60, 2.90, 2.24, 2.60, 1.42, 1.92, # Kr...
					   -7.26,  1.84, 1.12, 1.69, 0.94, 1.86, 1.66, # Tc... ???
						2.00,  1.00, 2.24, 1.51, 2.23, 1.12, 2.13, # Sn...
						1.22,  1.55, 0.71, 1.50, -7.96, 1.00, 0.51, # La...
						1.12, -0.10, 1.10, 0.26, 0.93, 0.00, 1.08, # Gd...
						0.76,  0.88, 0.13, 1.11, 0.27, 1.45, 1.35, # Lu...
						1.80,  1.01, 1.09, 0.90, 1.85, 0.71, -7.96, # Pt...
					   -7.96, -7.96, -7.96, -7.96, -7.96, 0.12, -7.96, # At...
					   -0.47, -7.96, -7.96, -7.96, -7.96, -7.96, -7.96, -7.96 # U...
					], dtype=np.float64)

atomweights = [
  {"H " :   1.008}, {"HE" :   4.003}, {"LI" :   6.939}, {"BE" :   9.013},
  {"B " :  10.810}, {"C " :  12.010}, {"N " :  14.010}, {"O " :  16.000},
  {"F " :  19.000}, {"NE" :  20.180}, {"NA" :  22.990}, {"MG" :  24.310},
  {"AL" :  26.980}, {"SI" :  28.090}, {"P " :  30.980}, {"S " :  32.070},
  {"CL" :  35.450}, {"AR" :  39.950}, {"K " :  39.100}, {"CA" :  40.080},
  {"SC" :  44.960}, {"TI" :  47.900}, {"V " :  50.940}, {"CR" :  52.000},
  {"MN" :  54.940}, {"FE" :  55.850}, {"CO" :  58.940}, {"NI" :  58.710},
  {"CU" :  63.550}, {"ZN" :  65.370}, {"GA" :  69.720}, {"GE" :  72.600},
  {"AS" :  74.920}, {"SE" :  78.960}, {"BR" :  79.910}, {"KR" :  83.800},
  {"RB" :  85.480}, {"SR" :  87.630}, {"Y " :  88.910}, {"ZR" :  91.220},
  {"NB" :  92.910}, {"MO" :  95.950}, {"TC" :  99.000}, {"RU" : 101.100},
  {"RH" : 102.900}, {"PD" : 106.400}, {"AG" : 107.900}, {"CD" : 112.400},
  {"IN" : 114.800}, {"SN" : 118.700}, {"SB" : 121.800}, {"TE" : 127.600},
  {"I " : 126.900}, {"XE" : 131.300}, {"CS" : 132.900}, {"BA" : 137.400},
  {"LA" : 138.900}, {"CE" : 140.100}, {"PR" : 140.900}, {"ND" : 144.300},
  {"PM" : 147.000}, {"SM" : 150.400}, {"EU" : 152.000}, {"GD" : 157.300},
  {"TB" : 158.900}, {"DY" : 162.500}, {"HO" : 164.900}, {"ER" : 167.300},
  {"TM" : 168.900}, {"YB" : 173.000}, {"LU" : 175.000}, {"HF" : 178.500},
  {"TA" : 181.000}, {"W " : 183.900}, {"RE" : 186.300}, {"OS" : 190.200},
  {"IR" : 192.200}, {"PT" : 195.100}, {"AU" : 197.000}, {"HG" : 200.600},
  {"TL" : 204.400}, {"PB" : 207.200}, {"BI" : 209.000}, {"PO" : 210.000},
  {"AT" : 211.000}, {"RN" : 222.000}, {"FR" : 223.000}, {"RA" : 226.100},
  {"AC" : 227.100}, {"TH" : 232.000}, {"PA" : 231.000}, {"U " : 238.000},
  {"NP" : 237.000}, {"PU" : 244.000}, {"AM" : 243.000}, {"CM" : 247.000},
  {"BK" : 247.000}, {"CF" : 251.000}, {"ES" : 254.000}
]

# @cython.auto_pickle(True)
cdef class Test:
	# from: https://stackoverflow.com/questions/36301322/pickle-cython-class-with-c-pointers
	cdef rh.Dusan me

	def __init__(self):
		self.me.n = 5
		self.me.dusan = 1.2

		self.me.run = <double *> malloc(self.me.n * sizeof(double))
		for idi in range(self.me.n):
			self.me.run[idi] = (idi+1)**2

		# self.me.again = rh.matrix_int(self.me.n, 10)
		# self.me.again = <int**> malloc(self.me.n * sizeof(int*))
		# for idi in range(self.me.n):
		# 	self.me.again[idi] = <int*> malloc(10 * sizeof(int))
		# self.me.again[0][0] = int(10)

	cdef bytes pack(self):
		a = <bytes>(<char*> self.me.run)[:sizeof(double)*self.me.n]
		# a = <bytes>(<char*> self.me.again)[:sizeof(int)*self.me.n*10]
		return a

	cdef void unpack(self, bytes data):
		tmp = <bytes>data[:sizeof(double)*self.me.n]
		memcpy(self.me.run, <char*>tmp, sizeof(double)*self.me.n)
		# tmp = <bytes>data#[sizeof(double)*self.me.n:]
		# memcpy(self.me.again, <char*>tmp, sizeof(int)*self.me.n*10)

	def get(self, idi):
		return self.me.run[idi]

	# def get_again(self):
	# 	return self.me.again[0][0]

	def __reduce__(self):
		data = self.pack()
		return (rebuild, (data,))

	# def __dealloc__(self):
	def __del__(self):
		free(self.me.run)
		# rh.freeMatrix(<void **>self.me.again)
		# if (self.me.again==NULL):
		# 	for idi in range(self.me.n):
		# 		free(self.me.again[idi])
		# 	free(self.me.again)

cpdef object rebuild(bytes data):
	c = Test()
	c.unpack(data)
	return c

cdef set_partition_functions(rh.Atmosphere *atmos):
	fp = open("../rh/Atoms/pf_Kurucz.input", "rb").read()
	buf = xdrlib.Unpacker(fp)

	cdef int Npf
	cdef cnp.ndarray[double, ndim=1, mode="c"] Tpf
	# cdef cnp.ndarray[double, ndim=2, mode="c"] pf
	cdef double[:, :] pf
	cdef cnp.ndarray[double, ndim=1, mode="c"] ionpot

	Npf = buf.unpack_int()
	Tpf = np.array(buf.unpack_farray(Npf, buf.unpack_double), dtype=np.float64)

	atmos.Npf = Npf
	atmos.Tpf = &Tpf[0]

	for ide in range(atmos.Nelem):
		pti = buf.unpack_int()
		Nstage = buf.unpack_int()
		atmos.elements[ide].Nstage = Nstage
		pf = np.array(buf.unpack_farray(Nstage*Npf, buf.unpack_double)).reshape(Nstage, Npf)
		pf = np.log(pf)
		ionpot = np.array(buf.unpack_farray(Nstage, buf.unpack_double))
		ionpot *= 19.8644746e-24

		atmos.elements[ide].pf = rh.matrix_double(Nstage, atmos.Npf)
		for ids in range(Nstage):
			atmos.elements[ide].pf[ids] = <double *>&pf[ids,0]

		atmos.elements[ide].ionpot = &ionpot[0]

cdef convert_2d(double **arr, Py_ssize_t nx, Py_ssize_t ny):
	cdef Py_ssize_t i
	cdef Py_ssize_t j
	cdef cnp.ndarray[cnp.float64_t, ndim=2] pyarr = np.zeros((nx, ny))
	for i in range(nx):
		for j in range(ny):
			pyarr[i,j] = arr[i][j]
	return pyarr

cpdef Pystring2char(lists):
	cdef int N = len(lists)
	cdef char* argv[140]

	py_string = [item.encode("utf-8") for item in lists]
	arr = (ctypes.c_char_p * N)(*py_string)
	for i_ in range(N):
		argv[i_] = arr[i_]

	return N, argv

class Populations(object):
	def __init__(self, ID, nlevel, nz, n, nstar):
		self.ID = ID
		self.nlevel = nlevel
		self.nz = nz
		self.n = n
		self.nstar = nstar

cdef rh.ne_solution get_solve_ne(_type):
	if _type=="NONE":
		return int(0)
	if _type=="ONCE":
		return int(1)
	if _type=="ITERATION":
		return int(2)
	
	sys.exit("Cannot assign 'solve_ne' attribute. Unsupported value.")

cdef rh.S_interpol_stokes get_S_interpol_stokes(_type):
	if _type=="DELO_PARABOLIC":
		return int(0)
	if _type=="DELO_BEZIER3":
		return int(1)

	sys.exit("Cannot assign 'S_interpolation_stokes' attribute. Unsupported value.")

cdef rh.S_interpol get_S_interpolation(_type):
	if _type=="S_LINEAR":
		return int(0)
	if _type=="S_PARABOLIC":
		return int(1)
	if _type=="S_BEZIER3":
		return int(2)

	sys.exit("Cannot assign 'S_interpolation' attribute. Unsupported value.")

cdef rh.StokesMode get_StokesMode(_type):
	if _type=="NO_STOKES":
		return int(0)
	if _type=="FIELD_FREE":
		return int(1)
	if _type=="POLARIZATION_FREE":
		return int(2)
	if _type=="FULL_STOKES":
		return int(3)

	sys.exit("Cannot assign 'StokesMode' attribute. Unsupported value.")

cdef rh.solution get_startJ(_type):
	if _type=="UNKNOWN":
		return int(-1)
	if _type=="LTE_POPULATIONS":
		return int(0)
	if _type=="ZERO_RADIATION":
		return int(1)
	if _type=="OLD_POPULATIONS":
		return int(2)
	if _type=="NEW_J":
		return int(3)
	if _type=="OLD_J":
		return int(4)

	sys.exit("Cannot assign 'startJ' attribute. Unsupported value.")

cdef int set_bool_value(flag):
	return int(1) if flag else int(0)

cdef int set_int_value(flag):
	return int(flag)

cdef class RH:
	cdef char* cwd[160]
	cdef rh.myRLK_Line lines
	cdef rh.InputData input
	cdef rh.Atmosphere atmos

	cdef public int Nwave
	cdef public double wavelength_min
	cdef public double wavelength_max
	cdef double* wavelength_vacuum

	def __init__(self):
		# convert the 'cwd' to the C char pointer
		py_list = cwd.split(" ")
		argc = len(py_list)
		py_string = [item.encode("utf-8") for item in py_list]
		arr = (ctypes.c_char_p * argc)(*py_string)
		for i_ in range(argc):
			self.cwd[i_] = arr[i_]

		self.atmos.Nelem = int(len(atomweights))
		
	def dummy(self):
		# just pass PF and Nstages trough the atmosphere and load it inside RH\
		# or create a separate function for this... I do not want to do this many times in inversion...
		rh.test_InputData(self.input, self.atmos)

	def set_abundances(self, abundances={}):
		self.input.abundances = <double *> malloc(self.atmos.Nelem * sizeof(double))
		for ida in range(self.atmos.Nelem):
			self.input.abundances[ida] = 10**(default_abundances[ida]-12)

		# update by the user provided values
		for key, value in abundances.items():
			self.input.abundances[key] = 10**(value-12)

		# multiply with metallicity and get total abundances
		self.atmos.totalAbund = self.input.abundances[0]
		self.atmos.wght_per_H = self.input.abundances[0] * list(atomweights[0].values())[0]
		for ida in range(1,self.atmos.Nelem):
			self.input.abundances[ida] *= self.input.metallicity
			self.atmos.totalAbund += self.input.abundances[ida] 
			self.atmos.wght_per_H += self.input.abundances[ida] * list(atomweights[ida].values())[0]

		self.atmos.avgMolWght = self.atmos.wght_per_H/self.atmos.totalAbund

	def set_keywords(self,
					 isum=-1, 
					 Ngdelay=0, 
					 Ngorder=0,
					 Ngperiod=1,
					 NmaxIter=1,
					 PRD_NmaxIter=1,
					 PRD_Ngdelay=0,
					 PRD_Ngorder=0,
					 PRD_Ngperiod=0,
					 NmaxScatter=0,
					 n_atomic_pars=0,
					 iterLimit=1e-2,
					 PRDiterLimit=1e-2,
					 metallicity=0.0,
					 solve_ne="NONE",
					 S_interpolation_stokes="DELO_BEZIER3",
					 S_interpolation="S_BEZIER3",
					 StokesMode="FULL_STOKES",
					 startJ="NEW_J",
					 magneto_optical=False,
					 PRD_angle_dep=False,
					 XRD=False,
					 Eddington=False,
					 backgr_pol=False,
					 allow_passive_bb=True,
					 NonICE=False,
					 rlkscatter=True,
					 xdr_endian=False,
					 old_background=False,
					 accelerate_mols=False,
					 do_fudge=False,
					 pyrhHSE=False,
					 get_atomic_rfs=False,
					 LS_Lande=True,
					 solve_NLTE=False,
					 verbose=False,
					 get_populations=False,
					 read_atom_model=False):
		#--- int attributes
		self.input.isum = set_int_value(isum)
		self.input.Ngdelay = set_int_value(Ngdelay)
		self.input.Ngorder = set_int_value(Ngorder)
		self.input.Ngperiod = set_int_value(Ngperiod)
		self.input.NmaxIter = set_int_value(NmaxIter)
		self.input.PRD_NmaxIter = set_int_value(PRD_NmaxIter)
		self.input.PRD_Ngdelay = set_int_value(PRD_Ngdelay)
		self.input.PRD_Ngorder = set_int_value(PRD_Ngorder)
		self.input.PRD_Ngperiod = set_int_value(PRD_Ngperiod)
		self.input.NmaxScatter = set_int_value(NmaxScatter)
		self.input.Nthreads = int(1) # always fixed
		self.input.n_atomic_pars = set_int_value(n_atomic_pars)

		#--- double attributes
		self.input.iterLimit = iterLimit
		self.input.PRDiterLimit = PRDiterLimit
		self.input.metallicity = 10**(metallicity)

		#--- enum attributes
		self.input.solve_ne = get_solve_ne(solve_ne)
		self.input.S_interpolation_stokes = get_S_interpol_stokes(S_interpolation_stokes)
		self.input.S_interpolation = get_S_interpolation(S_interpolation)
		self.input.StokesMode = get_StokesMode(StokesMode)
		self.input.startJ = get_startJ(startJ)

		#--- bool attributes
		self.input.magneto_optical = set_bool_value(magneto_optical)
		self.input.PRD_angle_dep = set_bool_value(PRD_angle_dep)
		self.input.XRD = set_bool_value(XRD)
		self.input.Eddington = set_bool_value(Eddington)
		self.input.backgr_pol = set_bool_value(backgr_pol)
		self.input.limit_memory = int(0)
		self.input.allow_passive_bb = set_bool_value(allow_passive_bb)
		self.input.NonICE = set_bool_value(NonICE)
		self.input.rlkscatter = set_bool_value(rlkscatter)
		self.input.xdr_endian = set_bool_value(xdr_endian)
		self.input.old_background = set_bool_value(old_background)
		self.input.accelerate_mols = set_bool_value(accelerate_mols)
		self.input.do_fudge = set_bool_value(do_fudge)
		self.input.pyrhHSE = set_bool_value(pyrhHSE)
		self.input.get_atomic_rfs = set_bool_value(get_atomic_rfs)
		self.input.LS_Lande = set_bool_value(LS_Lande)
		self.input.solve_NLTE = set_bool_value(solve_NLTE)
		self.input.verbose = set_bool_value(verbose)
		self.input.get_populations = set_bool_value(get_populations)
		self.input.read_atom_model = set_bool_value(read_atom_model)

	def set_elements(self):
		rh.set_elements(self.input, &self.atmos)
		set_partition_functions(&self.atmos)

	def get_RLK_lines(self):
		self.lines = rh.get_RLK_lines(self.cwd[0])

	def read_atom(self, atom_file_name):
		cdef char* pyrh_atom_file_name[100]
		
		py_list = atom_file_name.split(" ")
		argc = len(py_list)
		py_string = [item.encode("utf-8") for item in py_list]
		arr = (ctypes.c_char_p * argc)(*py_string)
		for i_ in range(argc):
			pyrh_atom_file_name[i_] = arr[i_]

		rh.read_atom_model(self.cwd[0], pyrh_atom_file_name[0])

@cython.boundscheck(False)
@cython.wraparound(False)
def get_ne_from_nH(cwd,
					 int atm_scale,
					 cnp.ndarray[double, ndim=1, mode="c"] scale,
					 cnp.ndarray[double, ndim=1, mode="c"] temperature,
					 cnp.ndarray[double, ndim=1, mode="c"] nH):
					 # cnp.ndarray[double, ndim=1, mode="c"] ne):
	cdef int Ndep = scale.size

	cdef char* argv[140]

	cdef cnp.ndarray[cnp.float64_t, ndim=1] ne = np.empty(Ndep)

	py_list = cwd.split(" ")
	argc = len(py_list)
	py_string = [item.encode("utf-8") for item in py_list]
	arr = (ctypes.c_char_p * argc)(*py_string)
	for i_ in range(argc):
		argv[i_] = arr[i_]

	rh.get_ne_from_nH(argv[0], 
			  atm_scale, 
			  Ndep, 
			  &scale[0], 
			  &temperature[0], 
			  &nH[0], 
			  &ne[0])

	return ne

@cython.boundscheck(False)
@cython.wraparound(False)
def hse(cwd,
		int atm_scale,
		cnp.ndarray[double, ndim=1, mode="c"] scale,
		cnp.ndarray[double, ndim=1, mode="c"] temp,
		double pg_top=0.1,
		cnp.ndarray[double, ndim=1, mode="c"] fudge_wave=None,
		cnp.ndarray[double, ndim=2, mode="c"] fudge_value=None,
		cnp.ndarray[int, ndim=1, mode="c"] atomic_number=None,
		cnp.ndarray[double, ndim=1, mode="c"] atomic_abundance=None,
		full_output=False):
	cdef int Ndep = scale.size

	cdef char* argv[140]

	cdef cnp.ndarray[cnp.float64_t, ndim=1] ne = np.empty(Ndep)
	cdef cnp.ndarray[cnp.float64_t, ndim=1] nHtot = np.empty(Ndep)
	cdef cnp.ndarray[cnp.float64_t, ndim=1] rho = np.empty(Ndep)
	cdef cnp.ndarray[cnp.float64_t, ndim=1] pg = np.empty(Ndep)
	pg[0] = pg_top

	py_list = cwd.split(" ")
	argc = len(py_list)
	py_string = [item.encode("utf-8") for item in py_list]
	arr = (ctypes.c_char_p * argc)(*py_string)
	for i_ in range(argc):
		argv[i_] = arr[i_]

	#--- fudge pointers
	cdef int fudge_num = 0
	cdef double* fudge_wave_ptr = NULL
	cdef double* fudge_value_ptr = NULL
	
	if (fudge_wave is not None) or (fudge_value is not None):
		fudge_num = fudge_wave.size
		fudge_wave_ptr = &fudge_wave[0]
		fudge_value_ptr = &fudge_value[0,0]

	#--- abundance pointers
	cdef int Nabun = 0
	cdef int* abundance_id_ptr = NULL
	cdef double* abundance_value_ptr = NULL

	if (atomic_number is not None) and (atomic_abundance is not None):
		Nabun = atomic_number.shape[0]
		if (Nabun!=atomic_abundance.shape[0]):
			print("\n  pyrh: Different number of 'atomic_number' and 'atomic_abundance'.\n")
			sys.exit()
		abundance_id_ptr = &atomic_number[0]
		abundance_value_ptr = &atomic_abundance[0]

	rh.hse(argv[0], Ndep,
			&scale[0], &temp[0],
			&ne[0], &nHtot[0], &rho[0], &pg[0],
			atm_scale,
			fudge_num, fudge_wave_ptr, fudge_value_ptr,
			Nabun, abundance_id_ptr, abundance_value_ptr)

	if full_output:
		return ne, nHtot, rho, pg

	del rho
	del pg

	return ne, nHtot

@cython.boundscheck(False)
@cython.wraparound(False)
def get_scales(cwd,
			  int atm_scale,
			  cnp.ndarray[double, ndim=1, mode="c"] scale,
			  cnp.ndarray[double, ndim=2, mode="c"] atmosphere,
			  double lam_ref,
			  cnp.ndarray[int, ndim=1, mode="c"] atomic_number=None,
			  cnp.ndarray[double, ndim=1, mode="c"] atomic_abundance=None):
	cdef int Ndep = atmosphere.shape[1]

	cdef char* argv[140]

	cdef cnp.ndarray[double, ndim=1, mode="c"] tau = np.empty(Ndep)
	cdef cnp.ndarray[double, ndim=1, mode="c"] height = np.empty(Ndep)
	cdef cnp.ndarray[double, ndim=1, mode="c"] cmass = np.empty(Ndep)

	py_list = cwd.split(" ")
	argc = len(py_list)
	py_string = [item.encode("utf-8") for item in py_list]
	arr = (ctypes.c_char_p * argc)(*py_string)
	for i_ in range(argc):
		argv[i_] = arr[i_]

	#--- abundance pointers
	cdef int Nabun = 0
	cdef int* abundance_id_ptr = NULL
	cdef double* abundance_value_ptr = NULL

	if (atomic_number is not None) and (atomic_abundance is not None):
		Nabun = atomic_number.shape[0]
		if (Nabun!=atomic_abundance.shape[0]):
			print("\n  pyrh: Different number of 'atomic_number' and 'atomic_abundance'.\n")
			sys.exit()
		abundance_id_ptr = &atomic_number[0]
		abundance_value_ptr = &atomic_abundance[0]

	rh.get_scales(argv[0], Ndep,
				 &scale[0], &atmosphere[1,0], 
				 &atmosphere[2,0], &atmosphere[3,0], &atmosphere[4,0],
				 &atmosphere[8,0], atm_scale,
				 lam_ref, &tau[0], &height[0], &cmass[0],
				 Nabun, abundance_id_ptr, abundance_value_ptr)

	return tau, height, cmass

@cython.boundscheck(False)
@cython.wraparound(False)
def compute1d(cwd,
				double mu,
				int atm_scale,
				cnp.ndarray[double, ndim=2, mode="c"] atmosphere,
				cnp.ndarray[double, ndim=1, mode="c"] wave,
				cnp.ndarray[int, ndim=1, mode="c"] loggf_ids=None,
				cnp.ndarray[double, ndim=1, mode="c"] loggf_values=None,
				cnp.ndarray[int, ndim=1, mode="c"] lam_ids=None,
				cnp.ndarray[double, ndim=1, mode="c"] lam_values=None,
				cnp.ndarray[double, ndim=1, mode="c"] fudge_wave=None,
				cnp.ndarray[double, ndim=2, mode="c"] fudge_value=None,
				cnp.ndarray[int, ndim=1, mode="c"] atomic_number=None,
				cnp.ndarray[double, ndim=1, mode="c"] atomic_abundance=None,
				get_atomic_rfs=False,
				get_populations=False):
	cdef int Ndep = atmosphere.shape[1]
	cdef int Nwave = wave.size
	
	#--- fudge pointers
	cdef int fudge_num = 0
	cdef double* fudge_wave_ptr = NULL
	cdef double* fudge_value_ptr = NULL
	
	if (fudge_wave is not None) or (fudge_value is not None):
		fudge_num = fudge_wave.size
		fudge_wave_ptr = &fudge_wave[0]
		fudge_value_ptr = &fudge_value[0,0]

	#--- log(gf) pointers
	cdef int Nloggf = 0
	cdef int* loggf_ids_ptr = NULL
	cdef double* loggf_values_ptr = NULL

	if (loggf_ids is not None) and (loggf_values is not None):
		Nloggf = loggf_ids.shape[0]
		if (Nloggf!=loggf_values.shape[0]):
			print("\n  pyrh: Different number of 'loggf_ids' and 'loggf_values'.\n")
			sys.exit()
		loggf_ids_ptr = &loggf_ids[0]
		loggf_values_ptr = &loggf_values[0]

	#--- lambda pointers
	cdef int Nlam = 0
	cdef int* lam_ids_ptr = NULL
	cdef double* lam_values_ptr = NULL

	if (lam_ids is not None) and (lam_values is not None):
		Nlam = lam_ids.shape[0]
		if (Nlam!=lam_values.shape[0]):
			print("\n  pyrh: Different number of 'lam_ids' and 'lam_values'.\n")
			sys.exit()
		lam_ids_ptr = &lam_ids[0]
		lam_values_ptr = &lam_values[0]

	#--- abundance pointers
	cdef int Nabun = 0
	cdef int* abundance_id_ptr = NULL
	cdef double* abundance_value_ptr = NULL

	if (atomic_number is not None) and (atomic_abundance is not None):
		Nabun = atomic_number.shape[0]
		if (Nabun!=atomic_abundance.shape[0]):
			print("\n  pyrh: Different number of 'atomic_number' and 'atomic_abundance'.\n")
			sys.exit()
		abundance_id_ptr = &atomic_number[0]
		abundance_value_ptr = &atomic_abundance[0]

	rh_get_atomic_rfs = 0
	if get_atomic_rfs:
		rh_get_atomic_rfs = 1
	
	rh_get_populations = 0
	if get_populations:
		rh_get_populations = 1
	
	cdef char* argv[140]

	py_list = cwd.split(" ")
	argc = len(py_list)
	py_string = [item.encode("utf-8") for item in py_list]
	arr = (ctypes.c_char_p * argc)(*py_string)
	for i_ in range(argc):
		argv[i_] = arr[i_]

	spec = rh.rhf1d(argv[0], mu, Ndep,
			 &atmosphere[0,0], &atmosphere[1,0], 
			 &atmosphere[2,0], &atmosphere[3,0], &atmosphere[4,0],
			 &atmosphere[5,0], &atmosphere[6,0], &atmosphere[7,0],
			 &atmosphere[8,0], atm_scale,
			 Nwave, &wave[0],
			 fudge_num, fudge_wave_ptr, fudge_value_ptr,
			 Nloggf, loggf_ids_ptr, loggf_values_ptr,
			 Nlam, lam_ids_ptr, lam_values_ptr,
			 Nabun, abundance_id_ptr, abundance_value_ptr,
			 rh_get_atomic_rfs, rh_get_populations,
			 0, argv[0])
			 # &self.wavetable[0], self.Nwave)

	lam = np.asarray(<cnp.float64_t[:spec.nlw]> spec.lam)
	sI = np.asarray(<cnp.float64_t[:spec.nlw]> spec.sI)

	sQ, sU, sV = None, None, None
	if spec.stokes:
		sQ = np.asarray(<cnp.float64_t[:spec.nlw]> spec.sQ)
		sU = np.asarray(<cnp.float64_t[:spec.nlw]> spec.sU)
		sV = np.asarray(<cnp.float64_t[:spec.nlw]> spec.sV)

	output = sI, sQ, sU, sV, lam

	if get_populations:
		populations = ()#[None]*spec.Nactive_atoms
		for ida in range(spec.Nactive_atoms):
			n = convert_2d(spec.atom_pops[ida].n, spec.atom_pops[ida].Nlevel, spec.atom_pops[ida].Nz)
			nstar = convert_2d(spec.atom_pops[ida].nstar, spec.atom_pops[ida].Nlevel, spec.atom_pops[ida].Nz)
			populations += Populations(
									ID=spec.atom_pops[ida].ID,
									nlevel=spec.atom_pops[ida].Nlevel,
									nz=spec.atom_pops[ida].Nz,
									n=n,
									nstar=nstar
									)
	if get_atomic_rfs:
		rf = convert_2d(spec.rfs, spec.nlw, Nloggf+Nlam)
		output = (output, rf.T)

	# to preserve the order of all output parameters
	if get_populations:
		output = (output, populations)

	# J = convert_2d(spec.J, spec.nlw, spec.Nrays)

	return output
