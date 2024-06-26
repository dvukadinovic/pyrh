Version 0.21:
-------------

-- updated the call signature for hse() method: fudge parameters are None by default.
-- added parameter 'full_output' to hse() method: by default returns only ne and nHtot, if 'full_output=True' return ne, nHtot, rho, pg.
-- updated rh/rh.h and rh/xdr.h to include header files from headers/types.h and headers/xdr.h instead from the system ones