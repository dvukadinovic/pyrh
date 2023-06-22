from setuptools import Extension, setup
from Cython.Build import cythonize
from distutils.sysconfig import get_config_vars as default_get_config_vars
import distutils.sysconfig as dsc
from platform import platform
import sys
import glob
import numpy

"""
Only CFLAGS needs to be modified
To Do:
  -- remove duplicate flags from CFLAGS (scarry, will skip it)
  -- if Mac OS X (with M1/M2 CPUs), remove the '-arch x86_64' flag from CFLAGS

"""

def remove_intel_arch(x):
    if type(x) is str:
        # if x=="-O2":
        #     return ""
        # if x.startswith("-O2 "):
        #     return remove_pthread(x[len("-O2 "):])
        # if x.endswith(" -O2"):
        #     return remove_pthread(x[:-len(" -O2")])
        return x.replace(" -arch x86_64 ", " ")
        # return x.replace(" -O3 ", " ")
    return x

def my_get_config_vars(*args):
    result = default_get_config_vars(*args)
    # sometimes result is a list and sometimes a dict
    if type(result) is list:
        return [remove_intel_arch(x) for x in result]
    elif type(result) is dict:
        return {k : remove_intel_arch(x) for k, x in result.items()}
    else:
        raise Exception("cannot handle type " + type(result))

if sys.platform=="macOS":
    dsc.get_config_vars = my_get_config_vars

print(dsc.get_config_vars()["CFLAGS"])

sys.exit()

rh_c_files = glob.glob("rh/*.c")

rhf1d = ["rh/rhf1d/anglequad.c", "rh/rhf1d/feautrier.c", "rh/rhf1d/multiatmos.c", \
	"rh/rhf1d/formal.c", "rh/rhf1d/piecestokes_1D.c", "rh/rhf1d/writeflux_xdr.c", \
	"rh/rhf1d/bezier_1D.c", "rh/rhf1d/hydrostat.c", "rh/rhf1d/piecewise_1D.c", "rh/rhf1d/riiplane.c",  \
	"rh/rhf1d/pyrh_compute1dray.c", "rh/rhf1d/pyrh_solveray.c",
	"rh/rhf1d/project.c", "rh/rhf1d/writegeom_xdr.c",
	"rh/rhf1d/pyrh_background.c", "rh/rhf1d/pyrh_hse.c"]

for item in rhf1d:
	rh_c_files.append(item)

rh_c_files.append("pyrh.pyx")
# rh_c_files.append("tools.pyx")

# for item in rh_c_files:
# 	print(item)

# remove some of the files
rh_c_files.remove("rh/collision_Oslo.c")
# rh_c_files.remove("rh/rhf1d/rf_ray.c")

setup(
	name="pyrh",
	version="0.2",
	author="Dusan Vukadionvic",
	ext_modules=cythonize([Extension("pyrh", rh_c_files,
							runtime_library_dirs=["rh"])],
							compiler_directives={"language_level" : "3"}),
	include_dirs=[numpy.get_include()]
)