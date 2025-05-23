Installation
============

Unix systems
------------

.. x86_64 arch systems
.. ^^^^^^^^^^^^^^^^^^^

To install the code you need to have ``Cython==0.29.33`` and ``numpy==1.24.2`` installed.

To install *pyrh*, simply run::

	python setup.py install

On every change of the source code, it needs to be recompiled since it is a shared library that needs to be updated are re-imported in the Python environment.

After installation, export the path to the ``pyrh/rh`` folder in the variable ``PYRH_PATH``. The best strategy is to set it directly in the ``~/.*rc`` file as::

	export PYRH_PATH="/absolte/path/to/the/pyrh/rh"

This path is used to locate atomic and molecular models (assumed to be in ``pyrh/rh/Atoms`` and ``pyrh/rh/Molecules``), abundance (in ``pyrh/rh/Atoms/abundance.input``) and the ABO coefficients (in ``pyrh/rh/Atoms/abo``).

.. caution::
	The instalation is tested also on a Mac PC with M1 chip. However, under the hood, RH uses SIMD instructions to compute a fast inversion of a 4x4 matrix, which is only available for x86_64 architectures. Therefore, expect 20-40% reduction in speed when using *pyrh* on Mac PCs with M1 chip (since it runs in arm64 architecture).

If everything passed without any errors, you have functional *pyrh* module! Now, lets compute some spectra!