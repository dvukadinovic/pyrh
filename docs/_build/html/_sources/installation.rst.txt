Installation
================================

For Unix based systems:

#. Install Cython module: https://cython.readthedocs.io/en/stable/src/quickstart/install.html. You would also need ``numpy`` module for `pyrh`.
#. Run ``python setup.py install`` and check for any possible errors. Make sure that OS is using the ``gcc`` compiler as default one. Otherwise you have to manually set the environment variable to point ot the correct ``gcc`` compiler.

If everything went without any error, you have functional `pyrh` module! Now, lets compute some spectra!