# Introduction

*pyrh* module is the Cython wrapper around [RH](https://github.com/han-uitenbroek/RH) synthesis code (Uitenbroek 2001). It allows access from Python environment to all RH methods to synthesise a spectrum in (non-)LTE, just like running RH code from a shell. The input parameters for RH rely still on standard `*.input` files. The user is suggested to consult the RH documentation for detailed overview of the code and usage of the input files.

# Installation

To install the code you need to have Cython==0.29.33 and numpy==1.24.2 installed.

To install *pyrh*, simply run:

	python setup.py install

On every change of the source code, it needs to be recompiled since it is a shared library that needs to be updated are re-imported in the Python environment.

After installation, export the path to the `pyrh/rh` folder in the variable `PYRH_PATH`. The best strategy is to set it directly in the `~/.*rc` file as:

	export PYRH_PATH="/absolte/path/to/the/pyrh/rh"

This path is used to locate atomic and molecular models (assumed to be in `pyrh/rh/Atoms` and `pyrh/rh/Molecules`), abundance (in `pyrh/rh/Atoms/abundance.input`) and the ABO coefficients (in `pyrh/rh/Atoms/abo`).

The instalation is tested also on Mac PC with M1 chip. However, under the hood, RH uses SIMD instructions to compute a fast inversion of a 4x4 matrix, which is only available for x86_64 architectures. Therefore, expect 20-40% reduction in speed when using *pyrh* on Mac PC with M1 chip (since it runs in arm64 architecture).

# Examples and test scripts

In directory `examples` you can find a simple script that shows how you can synthesise a spectrum using *pyrh*. 

In directory `tests` you can find more examples on how to call different methods available in *pyrh*. These are also used to verify changes made in the code and may not be user friendly.

# Important notes

We set the keyword `LIMIT_MEMORY` to be FALSE all the time. It is mannualy implemented inside the `pyrh_comute1dray.c` files. No threading is performed in wavelength domain as well.

# Issues

Please, report if there is any.