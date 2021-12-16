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