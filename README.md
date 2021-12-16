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
RLK lines ('atmos.Nrlk') we will enetr the function 'readKuruczLines()'. If we
already have them, no need to read because they will be copied into the memory 
where are they needed based on the input to 'rhf1d()' function.