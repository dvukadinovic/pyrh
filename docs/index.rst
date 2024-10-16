.. pyrh documentation master file, created by
   sphinx-quickstart on Mon May  1 19:42:29 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

*pyrh* documentation
====================

*pyrh* module is the Cython wrapper around [RH](https://github.com/han-uitenbroek/RH) NLTE synthesis code (Uitenbroek 2001). It allows access from Python environment to all RH methods to synthesise a spectrum in (non-)LTE, just like running RH code from a shell. The input parameters for RH rely still on standard `*.input` files. The user is suggested to consult the RH documentation for detailed overview of the code and usage of the input files.

Here, we will describe and show how to call the *pyrh* methods. Currently we can compute a spectrum in (non-)LTE, set an atmosphere in hydrostatic equilibrium (HSE) assuming an ideal gas equation and LTE populations of species, and compute the electron density from the known temperature and total hydrogen density.

*pyrh* is only called upon a single 1D atmospheric model, just like RH. For parallelised call you can use *globin* code which uses *pyrh*.

.. toctree::
   :maxdepth: 3
   :caption: Contents:

   installation
   quick_start
   input_data
   pyrh_api
.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`
