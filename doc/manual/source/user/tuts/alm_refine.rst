.. _tuts_alm_refine:

Wind-farm mesh refinement for Actuator Line simulation using Percept
====================================================================

This tutorial demonstrates the workflow for refining ABL meshes for use with
actuator line simulations using the Percept mesh adaptivity tool. We will start
with the precursor mesh and add nested zones of refinement around turbines of
interest so that the wakes are captured with adequate resolution necessary to
predict the impact on downstream turbine performance. We will perform the following steps

#. Use :program:`nalu_preprocess` to *tag* elements within the mesh that must be
   refined. In this exercise, we will perform two levels of refinement where the
   second level is nested within the first refinement zone. This step creates a
   ``turbine_refinement_field``, an element field, in the Exodus database. The
   refinement field is a scalar with a value ranging between 0 and 1. We will
   use this field as a threshold to control the regions where the refinement is
   applied by the :program:`mesh_adapt` utility in Percept.

#. Invoke Percept's :program:`mesh_adapt` utility twice to perform two levels of
   refinement. Each invocation will use the ``turbine_refinement_field``,
   created in the previous step, to determine the region where refinement is
   applied, the threshold is changed using YAML-formatted input files to
   :program:`mesh_adapt` during each call.


Prerequisites
-------------

To complete this tutorial you will need the Exodus mesh
(:file:`abl_1x1x1_10_mesh.exo`) generated in the :ref:`the previous tutorial
<tuts_abl_precursor>`. You will also need the input file for
:program:`nalu_preprocess` (:download:`abl_refine.yaml
<../files/tuts/abl_precursor/abl_refine.yaml>`)

Tag mesh regions for refinement
---------------------------------

In this step we will use :program:`nalu_preprocess` to create a refinement field
that will be used by :program:`mesh_adapt` to determine which elements are
selected for refinement. The input file that performs this action is shown below

.. literalinclude:: ../files/tuts/abl_precursor/abl_refine.yaml
   :lines: 14-30
   :linenos:

The mesh blocks targeted for refinement is provided as a list to the
``fluid_parts`` parameter (line 2), ``turbine_locations`` list the base
locations of the turbines in the wind farm that are being simulated,
``refinement_levels`` contain a list of length equal to the number of nested
refinement levels. Each entry in this list contains an array of four
non-dimensional lengths: the upstream, downstream, lateral, and vertical extent
of the refinement zones (as a multiple of rotor diameters) with respect to the
rotation center of the turbine. The orientation of the refinement boxes is
determined by the parameters provided within the ``orientation`` sub-dictionary.
In the current example, the boxes will be oriented along the wind direction
(:math:`245^\circ`) to match the ABL wind direction at hub-height used in the
previous tutorial.

.. note::

  It is recommended that the ``search_tolerance`` parameter in
  ``mesh_local_refinement`` section be set slightly larger than the coarset mesh
  resolution in the base ABL mesh chosen for refinement. This prevents jagged
  boundaries around the refinement zones as a result of roundoff and truncation
  errors. In our current example, this parameter was set to 11m based on the
  fact that the base mesh has a uniform resolution of 10m.


The output of :program:`nalu_preprocess` is shown below

::

  $ nalu_preprocess -i abl_refine.yaml

  Nalu Preprocessing Utility
  Input file: abl_refine.yaml
  Found 1 tasks
      - mesh_local_refinement

  Performing metadata updates...
  Metadata update completed
  Reading mesh bulk data... done.

  --------------------------------------------------
  Begin task: mesh_local_refinement
  Processing percept field: turbine_refinement_field
  Writing percept input files...
  	adapt1.yaml
  	adapt2.yaml
  Sample percept command line:
  mesh_adapt --refine=DEFAULT --input_mesh=mesh0.e --output_mesh=mesh1.e --RAR_info=adapt1.yaml
  End task: mesh_local_refinement

  All tasks completed; writing mesh...
  Exodus results file: mesh0.e

  Memory usage: Avg:  723.312 MB; Min:  723.312 MB; Max:  723.312 MB

Refine using Percept
--------------------

After executing :program:`nalu_preprocess` we should have :file:`mesh0.e`, the
Exodus database used as input for :program:`mesh_adapt` and two YAML files
:file:`adapt1.yaml` and :file:`adapt2.yaml` that contain the thresholds for each
level of refinement. To invoke Percept in serial mode, execute the following command

.. code-block:: bash

   # Refine the first level
   mesh_adapt --refine=DEFAULT --input_mesh=mesh0.e --output_mesh=mesh1.e --RAR_info=adapt1.yaml --progress_meter=1
   # Refine the second level
   mesh_adapt --refine=DEFAULT --input_mesh=mesh1.e --output_mesh=mesh2.e --RAR_info=adapt2.yaml --progress_meter=1

After successful execution of the two invocations of :program:`mesh_adapt`, the
refined mesh for use with actuator line wind farm simulations is saved in
:file:`mesh2.e`. Percept-based refinement creates pyramid and tetrahedral
elements at the refinement interfaces. These additional elements are added to
new mesh blocks (parts in STK parlance) that must be included in the Nalu input
file for simulation. Use :program:`ncdump` (see :ref:`previous tutorial
<tuts_abl_precursor>`) to examine the names of the new mesh blocks created by
Percept.

::

  $ ncdump -v eb_names mesh2.e
  #
  # OUTPUT TRUNCATED !!!
  #
  data:

   eb_names =
    "fluid",
    "fluid.pyramid_5._urpconv",
    "fluid.tetrahedron_4._urpconv",
    "fluid.pyramid_5._urpconv.Tetrahedron_4._urpconv" ;


For large meshes, parallel execution of Percept's :program:`mesh_adapt` utility
is recommended. A sample command line is shown below

.. code-block:: bash

   # Example mesh_adapt invocation in parallel.
   mpiexec -np 256 mesh_adapt --refine=DEFAULT --input_mesh=mesh0.e --output_mesh=mesh1.e --RAR_info=adapt2.yaml --progress_meter=1 --ioss_read_options="auto-decomp:yes" --ioss_write_options="large,auto-join:yes"

We pass ``auto-join:yes`` to IOSS write options so that the final mesh is
combined for subsequent use with a different number of MPI ranks with Nalu.

Troubleshooting tips
------------------------------

- Percept :program:`mesh_adapt` will hang if it runs out of memory without any
  error message. The user must ensure that enough memory is available to perform
  the refinements. Parallel execution on a larger number of nodes is the best
  solution to this problem.

- Percept creates long part names for the new mesh blocks it generates. These
  names are sometimes longer than the 32 characters allowed by SEACAS utilities
  for Exodus strings. Exodus mesh reading process will automatically truncate
  these names during read, but STK will throw an error if the full name is used
  to refer to the part. The user must take care to truncate the names to 32
  characters in the Nalu input file.

- Percept declares additional parts of form
  ``<BASE_PART>.pyramid_5._urpconv.Tetrahedron_4._urpconv`` in anticipation of
  possible refinement of pyramid elements into pyramids and tetrahedrons.
  However, the nested refinement strategy does not result in pyramids being
  refined and, therefore, this part remains empty. Currently, SEACAS and STK
  will throw an error if the user attempts to include this part in the Nalu
  input file during simulations.

- When using :program:`mesh_adapt` in parallel, appropriate IOSS read/write
  options must be specified to allow automatic decomposition of an undecomposed
  mesh and subsequent rejoin after parallel exection. Failure to provide
  appropriate options will lead to error during execution of
  :program:`mesh_adapt`.
