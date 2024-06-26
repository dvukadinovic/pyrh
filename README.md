# Installation

To install the package, simply run:

	python setup.py install

On every change of the code, it needs to be recompiled since it is a shared
library that needs to be changed are re-imported in the code in order to see
a change.

After installation, export the path to the `pyrh/rh` folder in variable
`PYRH_PATH`. The best strategy is to set it directly in the `~/.*rc` file as:

	export PYRH_PATH="/absolte/path/to/the/pyrh/rh"

This path is used to locate atomic and molecular models (assumed to be in
pyrh/rh/Atoms and pyrh/rh/Molecules), abundance (in pyrh/rh/Atoms) and for
ABO coefficients (in pyrh/rh/Atoms).

# Important notes

We set the keyword LIMIT_MEMORY to be FALSE all the time. It is mannualy
implemented inside the pyrh_comute1dray.c files.

# Changes

In `kurucz.c` in function `getUnsoldcross()`. The weight for Hydrogen atom 
was read from structure atmos, but now it is read from atmos.elements[0], which
point to the same object. This was necessary for independent reading of
Kurucz lines which will later be simply passed to the RH code to compute the
spectrum.

Added condition in `Background()` (from `background.c`). If there are no
RLK lines (`atmos.Nrlk=0`) we will enter the function `readKuruczLines()`. If we
already have them, no need to read because they will be copied into the memory 
where are they needed based on the input to `rhf1d()` function.

# Issues

Many...