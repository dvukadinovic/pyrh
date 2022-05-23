from setuptools import Extension, setup
from Cython.Build import cythonize
import glob
import numpy

rh_c_files = glob.glob("rh/*.c")
rh_c_files.remove("rh/collision_Oslo.c")
	

rhf1d = ["rh/rhf1d/anglequad.c", "rh/rhf1d/feautrier.c", "rh/rhf1d/multiatmos.c", \
	"rh/rhf1d/formal.c", "rh/rhf1d/piecestokes_1D.c", "rh/rhf1d/writeflux_xdr.c", \
	"rh/rhf1d/bezier_1D.c", "rh/rhf1d/hydrostat.c", "rh/rhf1d/piecewise_1D.c", "rh/rhf1d/riiplane.c",  \
	"rh/rhf1d/pyrh_compute1dray.c", "rh/rhf1d/pyrh_solveray.c",
	"rh/rhf1d/project.c", "rh/rhf1d/writegeom_xdr.c",
	"rh/rhf1d/pyrh_hse.c", "rh/rhf1d/pyrh_background.c"]

for item in rhf1d:
	rh_c_files.append(item)

rh_c_files.append("pyrh.pyx")
# rh_c_files.append("tools.pyx")

setup(
	name="pyrh",
	version="0.2",
	author="Dusan Vukadionvic",
	ext_modules=cythonize([Extension("pyrh", rh_c_files)],
							compiler_directives={'language_level' : "3"}),
	include_dirs=[numpy.get_include()]
)