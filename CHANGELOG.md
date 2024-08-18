Version 0.21:
-------------

-- update the call signature for hse() method: pg_top is now default parameter equal to 0.1 Pa.
-- updated the call signature for hse() method: fudge parameters are None by default.
-- updated the call signature for compute1d() method: added default parameters and adjusted the call to rh.rhf1d() method.
-- added parameter 'full_output' to hse() method: by default returns only ne and nHtot, if 'full_output=True' return ne, nHtot, rho, pg.
-- updated rh/rh.h and rh/xdr.h to include header files from headers/types.h and headers/xdr.h instead to look for the system ones.