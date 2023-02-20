# Installation

The RH code uses bool_t defined in rpc/types.h and is dependent on the
rpc/xdr.h IO (which is in pyrh obsolete; but is not removed from the code
completely). Before installation, types.h and xdr.h, located in ./headers, must be copied into /usr/include/rpc for RH to include them.

To install the package, do it as:

python setup.py install

On every change of the code, it needs to be recompiled since it is a shared
library that needs to be changed are re-imported in the code in order to see
change.

# Important notes

We set the keyword LIMIT_MEMORY to be FALSE all the time. It is mannualy
implemented inside the pyrh_comute1dray.c files.

# Changes to done

Added a keyword for a path to the Barklem table directory.

# Changes

In 'kurucz.c' in function 'getUnsoldcross()'. The weight for Hydrogen atom 
was read from structure atmos, but now it is read from atmos.elements[0], which
point to the same object. This was necessary for independent reading of
Kurucz lines which will later be simply passed to the RH code to compute the
spectrum.

Added condition in 'Background()' (from 'background.c'). If there are no
RLK lines ('atmos.Nrlk') we will enter the function 'readKuruczLines()'. If we
already have them, no need to read because they will be copied into the memory 
where are they needed based on the input to 'rhf1d()' function.

# Issues

argv[] list is not send to the rhf1d() in the format it should have as it is 
constructed in __init__().