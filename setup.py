from setuptools import Extension, setup
from Cython.Build import cythonize
import glob

rh_c_files = glob.glob("rh/*.c")
rh_c_files.remove("rh/collision_Oslo.c")
	

rhf1d = ["rh/rhf1d/anglequad.c", "rh/rhf1d/feautrier.c", "rh/rhf1d/multiatmos.c", \
	"rh/rhf1d/formal.c", "rh/rhf1d/piecestokes_1D.c", "rh/rhf1d/writeflux_xdr.c", \
	"rh/rhf1d/bezier_1D.c", "rh/rhf1d/hydrostat.c", "rh/rhf1d/piecewise_1D.c", "rh/rhf1d/riiplane.c",  \
	"rh/rhf1d/dv_rhf1d.c", "rh/rhf1d/project.c"]

"""
had to delete:
rh/rhf1d/writegeom_xdr.c --> multiple topology definitions

"""

for item in rhf1d:
	rh_c_files.append(item)

rh_c_files.append("pyrh.pyx")
# rh_c_files.append("tools.pyx")

# rh_c_files.append("rh/options.c")
# rh_c_files.append("rh/parse.c")

setup(
	name="pyrh",
	version="0.1",
	author="Dusan Vukadionvic",
	ext_modules=cythonize([Extension("pyrh", rh_c_files)])
)