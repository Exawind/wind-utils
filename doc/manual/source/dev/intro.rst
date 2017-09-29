.. _dev_intro:

Introduction
============

This part of the documentation is intended for users who wish to extend or add
new functionality to the NaluWindUtilities toolsuite. End users who want to use
existing utilities should consult the :ref:`user_guide` for documentation on
standalone utilities.

Version Control System
----------------------

Like `Nalu <http://nalu.readthedocs.io/en/latest/>`_, NaluWindUtils uses `Git
SCM <https://www.git-scm.com>`_ to track all development activity. All
development is coordinated through the `Github repository
<https://github.com/NaluCFD/NaluWindUtils>`_. `Pro Git
<https://www.git-scm.com/book/en/v2>`_, a book that covers all aspects of Git is
a good resource for users unfamiliar with Git SCM. `Github Desktop
<https://desktop.github.com>`_ and `Git Kraken <https://www.gitkraken.com>`_ are
two options for users who prefer a GUI based interaction with Git source code.

.. _dev_docs_build:

Building API Documentation
--------------------------

In-source comments can be compiled and viewed as HTML files using `Doxygen
<http://www.stack.nl/~dimitri/doxygen/index.html>`_. If you want to generate
class inheritance and other collaboration diagrams, then you will need to
install `Graphviz <http://www.graphviz.org>`_ in addition to Doxygen.

#. API Documentation generation is disabled by default in CMake. Users will have
   to enable this by turning on the :cmakeval:`ENABLE_DOXYGEN_DOCS` flag.

#. Run ``make api-docs`` to generate the documentation in HTML form.

The resulting documentation will be available in :file:`doc/doxygen/html/`
within the CMake build directory.

Contributing
-------------

The project welcomes contributions from the wind research community. Users can
contribute to the source code using the normal `Github fork and pull request
workflow <https://guides.github.com/activities/forking/>`_. Please follow these
general guidelines when submitting pull requests to this project

* All C++ code must conform to the C++11 standard. Consult `C++ Core Guidelines
  <http://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines>`_ on
  best-practices to writing idiomatic C++ code.

* Check and fix all compiler warnings before submitting pull requests. Use
  ``-Wall -Wextra -pedantic`` options with GNU GCC or LLVM/Clang to check for
  warnings.

* New feature pull-requests must include doxygen-compatible *in source*
  documentation, additions to user manual describing the enchancements and their
  usage, as well as the necessary updates to CMake files to enable configuration
  and build of these capabilities.

* Prefer Markdown format when documenting code using Doxgen-compatible comments.

* Avoid incurring additional third-party library (TPL) dependencies beyond what
  is required for building Nalu. In cases where this is unavoidable, please
  discuss this with the development team by creating an issue on `issues page
  <https://github.com/NaluCFD/NaluWindUtils/issues>`_ before submitting the pull
  request.
