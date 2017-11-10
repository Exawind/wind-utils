.. _util_nalu_postprocess_exe:

``nalu_postprocess`` -- Nalu Post-processing Utilities
======================================================

This utility loads an Exodus-II solution file and performs various
post-processing *tasks* on the database. Currently, the following *tasks* have
been implemented within this utility.

==================== ===================================================
Task type            Description
==================== ===================================================
``abl_statistics``   Calculate various ABL statistics of interest
==================== ===================================================

The input file (:download:`download <files/nalu_postprocess.yaml>`) must contain
a **nalu_postprocess** section a shown below. Input options for various *tasks*
are provided as sub-sections within **nalu_postprocess** with the corresponding
task names under tasks.

.. literalinclude:: files/nalu_postprocess.yaml
   :language: yaml

Command line invocation
-----------------------

.. code-block:: bash

   mpiexec -np <N> nalu_postprocess -i [YAML_INPUT_FILE]

.. program:: nalu_postprocess

.. option:: -i, --input-file

   Name of the YAML input file to be used. Default: ``nalu_postprocess.yaml``.

Common input file options
-------------------------

.. confval:: input_db

   Path to an existing Exodus-II mesh database file, e.g.., ``ablPrecursor.e``

.. confval:: tasks

   A list of task names that define the various pre-processing tasks that will
   be performed on the input mesh database by this utility. The program expects
   to find additional sections with headings matching the task names that
   provide additional inputs for individual tasks.


``abl_statistics``
------------------

This task computes various various statistics relevant for ABL simulations and
outputs vertical profiles of various quantities of interest.

.. literalinclude:: files/nalu_postprocess.yaml
   :language: yaml
   :lines: 12-25
