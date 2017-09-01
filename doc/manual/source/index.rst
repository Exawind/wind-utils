
###############################
Nalu Wind Utilities User Manual
###############################

.. only:: html

   :Version: |version|
   :Date: |today|

NaluWindUtils is a companion software library to `Nalu
<http://nalu.readthedocs.io/en/latest/>`_ --- a generalized, unstructured,
massively parallel, low-Mach flow solver for wind energy applications. As the
name indicates, this software repository provides various meshing, pre- and
post-processing utilities for use with the Nalu CFD code to aid setup and
analysis of wind energy LES problems. This software is licensed under `Apache
License Version 2.0 <http://www.apache.org/licenses/LICENSE-2.0>`_ open-source
license.

The source code is hosted and all development is coordinated through the `Github
repository <https://github.com/NaluCFD/NaluWindUtils>`_ under the `NaluCFD
organization <https://github.com/NaluCFD>`_ umbrella. The official documentation
for all released and development versions are hosted on `ReadTheDocs
<http://naluwindutils.readthedocs.io/en/latest/>`_. Users are welcome to submit
issues, bugs, or questions via the `issues page
<https://github.com/NaluCFD/NaluWindUtils/issues>`_. Users are also encouraged
to contribute to the source code and documentation using `pull requests
<https://github.com/NaluCFD/NaluWindUtils/pulls>`_ using the normal `Github fork
and pull request workflow <https://guides.github.com/activities/forking/>`_.

This documentation is divided into two parts:

:ref:`user_guide`

   Directed towards end-users, this part provides detailed documentation
   regarding installation and usage of the various utilities available within
   this library. Here you will find a comprehensive listing of all available
   utilties, and information regarding their usage and current
   limitations that the users must be aware of.

:ref:`dev_guide`

   The developer guide is targeted towards users wishing to extend the
   functionality provided within this library. Here you will find details
   regarding the code structure, API supported by various classes, and links to
   source code documentation extracted using Doxygen.

**Acknowledgements**

This software is developed by researchers at `NREL <https://www.nrel.gov>`_ and
`Sandia National Laboratories <http://www.sandia.gov>`_ with funding from DOE's
`Exascale Computing Project <https://exascaleproject.org>`_ and DOE WETO
`Atmosphere to electrons (A2e) <https://a2e.energy.gov>`_ research initiative.

.. toctree::
   :maxdepth: 3

   user/index
   dev/index

