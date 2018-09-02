.. _tuts_abl_precursor:

Pre-processing for ABL precursor runs
=====================================

This tutorial walks through the steps required to create an ABL mesh and
initialize the fields for an ABL precursor run. In this tutorial, you will use
the :ref:`abl_mesh <util_abl_mesh_exe>` and certain capabilities of
:ref:`nalu_preprocess <util_nalu_preprocess_exe>`. The steps covered in this tutorial are

   #. Generate a :math:`1 \times 1 \times 1` km HEX block mesh with uniform
      resolution of 10m in all three directions.

   #. Generate a *sampling plane* at hub-height (90m) where the velocity field
      will be sampled to force it to a desired wind speed and direction using a
      driving pressure gradient source term.

   #. Initialize the velocity and temperature field to the desired profile as a
      function of height, and add perturbations to the fields to kick-off
      turbulence generation.


Prerequisites
--------------

Before attempting this tutorial, you should have a compiled version of
NaluWindUtils. Please consult the :ref:`user_installation` section to fetch,
configure, and compile the latest version of the source code. You can also
download the input file (:download:`abl_setup.yaml
<../files/tuts/abl_precursor/abl_setup.yaml>`) that will be used with
:program:`abl_mesh` and :program:`nalu_preprocess` executables.


Generate ABL precursor mesh
---------------------------

In this step, we will use the :program:`abl_mesh` utility to generate :math:`1
\times 1 \times 1` km with a uniform resolution of 10m in all three directions.
The domain will span :math:`[0, 1000]` m in each direction. The relevant section
in the input file is shown below

.. literalinclude:: ../files/tuts/abl_precursor/abl_setup.yaml
   :lines: 6-20
   :linenos:

With this section saved in the input file :file:`abl_setup.yaml`, the sample
interaction is shown below

::

  $ abl_mesh -i abl_setup.yaml

  Nalu ABL Mesh Generation Utility
  Input file: abl_setup.yaml
  HexBlockBase: Registering parts to meta data
  	Mesh block: fluid
  Num. nodes = 1030301; Num elements = 1000000
  	Generating node IDs...
  	Creating nodes... 10% 20% 30% 40% 50% 60% 70% 80% 90%
  	Generating element IDs...
  	Creating elements... 10% 20% 30% 40% 50% 60% 70% 80% 90%
  	Finalizing bulk modifications...
  	Generating X Sideset: west
  	Generating X Sideset: east
  	Generating Y Sideset: south
  	Generating Y Sideset: north
  	Generating Z Sideset: terrain
  	Generating Z Sideset: top
  	Generating coordinates...
  	 Generating x spacing: constant_spacing
  	 Generating y spacing: constant_spacing
  	 Generating z spacing: constant_spacing
  Writing mesh to file: abl_1x1x1_10_mesh.exo

  Memory usage: Avg:  553.148 MB; Min:  553.148 MB; Max:  553.148 MB


Initializing fields and sampling planes
---------------------------------------

In the next step we will use :program:`nalu_preprocess` to setup the fields
necessary for a precursor simulation. The relevant section of the input file is
shown below

.. literalinclude:: ../files/tuts/abl_precursor/abl_setup.yaml
   :lines: 21-58
   :linenos:

The following actions are performed

#. Lines 14--18: Initialize a constant velocity field such that the wind speed is 8.0 m/s
   along :math:`245^\circ` compass direction.

#. Lines 24--26: A constant temperature field of 300K till 650m and then a
   capping inversion between 650m to 750m and a temperature gradient of 0.003
   K/m above the capping inversion zone.

#. Pertubations to the velocity (lines 19--22) and temperature field (lines
   27--30) to kick off turbulence generation during the precursor run. The
   velocity field perturbations are similar to those generated in SOWFA for ABL
   precursor runs.

#. The mesh generated in the previous step is used as input (line 5), and a new
   file is written out with the new fields and the sampling plane (line 6).

Output from the execution of :program:`nalu_preprocess` with this input file is shown below

::

  $ nalu_preprocess -i abl_setup.yaml

  Nalu Preprocessing Utility
  Input file: abl_setup.yaml
  Found 1 tasks
      - init_abl_fields

  Performing metadata updates...
  Metadata update completed
  Reading mesh bulk data... done.

  --------------------------------------------------
  Begin task: init_abl_fields
  Generating ABL fields
  End task: init_abl_fields

  All tasks completed; writing mesh...
  Exodus results file: abl_1x1x1_10.exo

  Memory usage: Avg:  786.082 MB; Min:  786.082 MB; Max:  786.082 MB

Using ``ncdump`` to examine mesh metadata
-----------------------------------------

:program:`ncdump` is a NetCDF utility that is built and installed as a depedency
of Trilinos. Since Trilinos is a dependency of NaluWindUtils, you should have
:program:`ncdump` available in your path if Trilinos and its dependencies were
loaded properly (either via ``spack`` or ``module load``). :program:`ncdump` is
useful to quickly examine the Exodus file metadata from the command line. Invoke
the command with ``-h`` option to quickly see the number of nodes and elements
in a mesh

::

   $ ncdump -h abl_1x1x1_10.exo
   netcdf abl_1x1x1_10 {
   dimensions:
       len_string = 33 ;
       len_line = 81 ;
       four = 4 ;
       num_qa_rec = 1 ;
       num_info = 2 ;
       len_name = 33 ;
       num_dim = 3 ;
       time_step = UNLIMITED ; // (1 currently)
       num_nodes = 1040502 ;
       num_elem = 1000000 ;
       num_el_blk = 1 ;
       num_node_sets = 1 ;
       num_side_sets = 6 ;
       num_el_in_blk1 = 1000000 ;
       num_nod_per_el1 = 8 ;
       num_nod_ns1 = 10201 ;
       num_side_ss1 = 10000 ;
       num_df_ss1 = 40000 ;
       num_side_ss2 = 10000 ;
       num_df_ss2 = 40000 ;
       num_side_ss3 = 10000 ;
       num_df_ss3 = 40000 ;
       num_side_ss4 = 10000 ;
       num_df_ss4 = 40000 ;
       num_side_ss5 = 10000 ;
       num_df_ss5 = 40000 ;
       num_side_ss6 = 10000 ;
       num_df_ss6 = 40000 ;
       num_nod_var = 4 ;

For the ABL precursor mesh generated using :program:`abl_mesh` we have 1 mesh
block (``num_el_blk``) that has one million elements (``num_el_in_blk1``)
composed of Hexahedral elements with 8 nodes per element (``num_nod_per_el1``).
There are 4 nodal field variables (``num_nod_var``) stored in this database that
were created by :program:`nalu_preprocess` utility. Finally, there are 6
sidesets (``num_side_sets``) each with 10,000 faces, and one node set
(``num_node_sets``) that contains 10201 nodes that were created as a sampling
plane at hub height of 90m during the pre-processing step.

Use the ``-v`` flag with the desired variable names (separated by commas) to
examine the contents of those variables. For example, to output the mesh blocks
(``eb_names``), sidesets or boundaries (``ss_names``), nodal (``name_nod_var``)
and element fields (``name_elem_var``) present in an Exodus database:

::

  $ ncdump -v eb_names,ss_names,name_nod_var abl_1x1x1_10.exo
  #
  # OUTPUT TRUNCATED !!!
  #
  data:

   eb_names =
    "fluid" ;

   ss_names =
    "west",
    "east",
    "south",
    "north",
    "terrain",
    "top" ;

   name_nod_var =
    "temperature",
    "velocity_x",
    "velocity_y",
    "velocity_z" ;

As seen in the output for ``name_nod_var``, Exodus file contains
``temperature``, a scalar field, and ``velocity``, a vector field. Internally,
exodus stores each component of a vector or tensor field as a separate variable.
The mesh block is called ``fluid`` can should be referred as such in the
pre-processing tasks or within the Nalu input file. As indicated in the
``dimensions``, this file contains one ``time_step``, you can use ``-v
time_whole`` to determine the timesteps that are currently stored in the Exodus
database.
