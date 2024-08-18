# Introduction

*pyrh* module is the Cython wrapper around RH synthesis code (Uitenbroek 2001; https://github.com/han-uitenbroek/RH). It allows access from Python environment to all RH methods to synthesise a spectrum in (non-)LTE, just like running RH code from a shell. The input parameters for RH rely still on standard `*.input` files. The user is suggested to consult the RH documentation for detailed overview of the code and usage of the input files.

# Installation

To install the code you need to have Cython==0.29.33 and numpy==1.24.2 installed.

To install *pyrh*, simply run:

	python setup.py install

On every change of the code, it needs to be recompiled since it is a shared library that needs to be updated are re-imported in the code in order to see a change made in the source code

After installation, export the path to the `pyrh/rh` folder in variable `PYRH_PATH`. The best strategy is to set it directly in the `~/.*rc` file as:

	export PYRH_PATH="/absolte/path/to/the/pyrh/rh"

This path is used to locate atomic and molecular models (assumed to be in `pyrh/rh/Atoms` and `pyrh/rh/Molecules`), abundance (in `pyrh/rh/Atoms/abundance.input`) and for the ABO coefficients (in `pyrh/rh/Atoms/abo`).

# Examples and test scripts

In directory `examples` you can find a simple script that shows how you can synthesise a spectrum using *pyrh*. 

In directory `tests` you can find more examples on how to call different methods available in *pyrh*. These are also used to verify changes made in the code and may not be user friendly.

# Important notes

We set the keyword LIMIT_MEMORY to be FALSE all the time. It is mannualy implemented inside the pyrh_comute1dray.c files. No threading is performed in wavelength domain.

# Issues

Many...