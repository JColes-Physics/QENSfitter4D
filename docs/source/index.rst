.. QENSfitter4D documentation master file, created by
   sphinx-quickstart on Fri Jul 24 18:11:52 2026.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

QENSfitter4D documentation
==========================

QENSfitter4D is a NeXpy plugin which provides a GUI resource for analyzing 4D-QENS data in the nexus format.

This package requires the use of the three following packages It is usefull to familiarize oneself with at least the first two packages, NeXpy and NeXus as these will be part of the user interactions. LMFIT is usefull for understanding how the fitting algorithms achieve their results.

**NeXpy**
  `NeXpy <https://github.com/nexpy/nexpy>`__ provides the GUI
  interface for loading, inspecting, plotting, and manipulating NeXus
  data, with an embedded IPython shell and script editor.

  .. image:: https://img.shields.io/pypi/v/nexpy.svg
     :target: https://pypi.python.org/pypi/nexpy

  .. image:: https://img.shields.io/conda/vn/conda-forge/nexpy
     :target: https://anaconda.org/conda-forge/nexpy

**nexusformat**
  The API for reading, modifying, and writing NeXus data is provided by
  the `nexusformat <https://github.com/nexpy/nexusformat>`__ package,
  which utilizes `h5py <http://www.h5py.org/>`__ for loading and saving
  the data in HDF5 files.

  .. image:: https://img.shields.io/pypi/v/nexusformat.svg
     :target: https://pypi.python.org/pypi/nexusformat

  .. image:: https://img.shields.io/conda/vn/conda-forge/nexusformat
     :target: https://anaconda.org/conda-forge/nexusformat

**LMFIT**
   The API for performing automated fitting functions for dynamic structure factor is
   provided by the `LMFit <https://lmfit.github.io/lmfit-py/>`__ package.

   .. image:: https://img.shields.io/pypi/v/lmfit
      :target: https://pypi.org/project/lmfit/

   .. image:: https://img.shields.io/conda/vn/conda-forge/lmfit
      :target: https://anaconda.org/channels/conda-forge/packages/lmfit/overview


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   includeme
   fileformatting
   initializing_input_params
   rasterscanning
   second_voigt

.. toctree::
   :maxdepth: 2

   api

