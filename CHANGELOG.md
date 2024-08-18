Version 0.21:
-------------

* updated the call signature for `hse()` method: pg_top is now default parameter equal to 0.1 Pa.
* updated the call signature for `hse()` method: fudge parameters are None by default.
* updated the call signature for `compute1d()` method: added default parameters and adjusted the call to `rh.rhf1d()` method.
* added parameter `full_output` to `hse()` method: by default returns only `ne` and `nHtot`, if `full_output=True` return `ne, nHtot, rho, pg`.
* updated `rh/rh.h` and `rh/xdr.h` to include header files from `headers/types.h` and `headers/xdr.h` instead to look for the system ones.
* reverted change for conditional reading a Kurucz line list (from version 0.10).

Version 0.10:
-------------

* Added condition in `Background()` (from `background.c`). If there are no RLK lines (`atmos.Nrlk=0`) we will enter the function `readKuruczLines()`. If we already have them, no need to read because they will be copied into the memory where are they needed based on the input to `rhf1d()` function.
* In `kurucz.c` in function `getUnsoldcross()`. The weight for hydrogen atom was read from structure atmos, but now it is read from `atmos.elements[0]`, which point to the same object. This was necessary for independent reading of Kurucz lines which will later be simply passed to the RH code to compute the spectrum.
