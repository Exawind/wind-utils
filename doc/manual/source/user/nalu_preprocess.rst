.. _util_nalu_preprocess_exe:


``nalu_preprocess`` -- Nalu Preprocessing Utilities
===================================================

This utility loads an input mesh and performs various pre-processing *tasks* so
that the resulting output database can be used in a wind LES simulation.
Currently, the following *tasks* have been implemented within this utility.

=========================  ===========================================================
Task type                  Description
=========================  ===========================================================
``init_abl_fields``        Initialize ABL velocity and temperature fields
``init_channel_fields``    Initialize channel velocity fields
``generate_planes``        Generate horizontal sampling planes for ``dp/dx`` forcing
``create_bdy_io_mesh``     Create an I/O transfer mesh for sampling inflow planes
``mesh_local_refinement``  Local refinement around turbines for wind farm simulations
``rotate_mesh``            Rotate mesh
=========================  ===========================================================

.. warning::

   Not all tasks are capable of running in parallel. Please consult
   documentation of individual tasks to determine if it is safe to run it in
   parallel using MPI. It might be necessary to set
   :confval:`automatic_decomposition_type` when running in parallel.

The input file (:download:`download <files/nalu_preprocess.yaml>`) must contain
a **nalu_preprocess** section as shown below. Input options for the individual
tasks are provided as sub-sections within **nalu_preprocess** with the
corresponding task names provided under ``tasks``. For example, in the sample
shown below, the program will expect to see two sub-sections, namely
``init_abl_fields`` and ``generate_planes`` based on the list of tasks shown in
lines 22-23.

.. literalinclude:: files/nalu_preprocess.yaml
   :language: yaml
   :linenos:
   :lines: 1-23

Command line invocation
-----------------------

.. code-block:: bash

   mpiexec -np <N> nalu_preprocess -i [YAML_INPUT_FILE]

.. program:: nalu_preprocess

.. option:: -i, --input-file

   Name of the YAML input file to be used. Default: ``nalu_preprocess.yaml``.

Common input file options
-------------------------

.. confval:: input_db

   Path to an existing Exodus-II mesh database file, e.g.., ``ablNeutralMesh.g``

.. confval:: output_db

   Filename where the pre-processed results database is output, e.g.,
   ``ablNeutralPrecursor.g``

.. confval:: automatic_decomposition_type

   Used only for parallel runs, this indicates how the a single mesh database
   must be decomposed amongst the MPI processes during initialization. This
   option should not be used if the mesh has already been decomposed by an
   external utility. Possible values are:

   ==========  ==========================================================
   Value       Description
   ==========  ==========================================================
   rcb         recursive coordinate bisection
   rib         recursive inertial bisection
   linear      elements in order first n/p to proc 0, next to proc 1.
   cyclic      elements handed out to id % proc_count
   ==========  ==========================================================

.. confval:: tasks

   A list of task names that define the various pre-processing tasks that will
   be performed on the input mesh database by this utility. The program expects
   to find additional sections with headings matching the task names that
   provide additional inputs for individual tasks. By default, the task names
   found within the list should correspond to one of the **task types**
   discussed earlier in this section. If the user desires to use custom names,
   then the exact task type should be provided with a ``type`` within the task
   section. A specific use-case where this is useful is when the user desires to
   rotate the mesh, perform additional operations, and, finally, rotate it back
   to the original orientation.

   .. code-block:: yaml
      :linenos:
      :emphasize-lines: 2,4,6,7,15,16

      tasks:
        - rotate_mesh_ccw  # Rotate mesh such that sides align with XYZ axes
        - generate_planes  # Generate sampling planes using bounding box
        - rotate_mesh_cw   # Rotate mesh back to the original orientation

      rotate_mesh_ccw:
        task_type: rotate_mesh
        mesh_parts:
          - unspecified-2-hex

        angle: 30.0
        origin: [500.0, 0.0, 0.0]
        axis: [0.0, 0.0, 1.0]

      rotate_mesh_cw:
        task_type: rotate_mesh
        mesh_parts:
          - unspecified-2-hex
          - zplane_0080.0         # Rotate auto generated parts also

        angle: -30.0
        origin: [500.0, 0.0, 0.0]
        axis: [0.0, 0.0, 1.0]

.. confval:: transfer_fields

   A Boolean flag indicating whether the time histories of the fields available
   in the input mesh database must be transferred to the output database.
   Default: ``false``.

.. confval:: ioss_8bit_ints

   A Boolean flag indicating whether the output database must be written out
   with 8-bit integer support. Default: ``false``.

``init_abl_fields``
-------------------

This task initializes the vertical velocity and temperature profiles for use
with an ABL precursor simulations based on the parameters provided by the user
and writes it out to the :confval:`output_db`. It is safe to run
``init_abl_fields`` in parallel. A sample invocation is shown below

.. code-block:: yaml
   :linenos:

   init_abl_fields:
     fluid_parts: [fluid]

     temperature:
       heights: [    0, 650.0, 750.0, 10750.0]
       values:  [280.0, 280.0, 288.0,   318.0]

       # Optional section to add random perturbations to temperature field
       perturbations:
         amplitude: 0.8 # in Kelvin
         cutoff_height: 600.0 # Perturbations below capping inversion
         skip_periodic_parts: [east, west, north, south]

     velocity:
       heights: [0.0, 10.0, 30.0, 70.0, 100.0, 650.0, 10000.0]
       values:
         - [ 0.0, 0.0, 0.0]
         - [4.81947, -4.81947, 0.0]
         - [5.63845, -5.63845, 0.0]
         - [6.36396, -6.36396, 0.0]
         - [6.69663, -6.69663, 0.0]
         - [8.74957, -8.74957, 0.0]
         - [8.74957, -8.74957, 0.0]

       # Optional section to add sinusoidal streaks to the velocity field
       perturbations:
         reference_height: 50.0   # Reference height for damping
         amplitude: [1.0, 1.0]    # Perturbation amplitudes in Ux and Uy
         periods: [4.0, 4.0]      # Num. periods in x and y directions

.. confval:: fluid_parts

   A list of element block names where the velocity and/or temperature fields
   are to be initialized.

.. confval:: temperature

   A YAML dictionary containing two arrays: ``heights`` and the corresponding
   ``values`` at those heights. The data must be provided in SI units. No
   conversion is performed within the code.

   The temperature section can contain an optional section ``perturbations``
   (lines 8-12) that will add fluctuations to the temperature field. It requires
   three parameters: 1. the amplitude of oscillations (in degrees Kelvin), 2.
   the cutoff height above which perturbations are not added, and a list of
   sidesets where the perturbations should not be added. It is important that
   the perturbations are not added to the periodic sidesets, otherwise the Nalu
   simulations will show spurious flow structures.

.. confval:: velocity

   A YAML dictionary containing two arrays: ``heights`` and the corresponding
   ``values`` at those heights. The data must be provided in SI units. No
   conversion is performed within the code. The values in this case are two
   dimensional lists of shape ``[nheights, 3]`` where ``nheights`` is the length
   of the `heights` array provided.

   Like temperature, the user can add sinusoidal streaks to the velocity field
   to trigger the turbulence generation -- see lines 25-29. The implementation
   follows the method used in SOWFA.

.. note::

   Only one of the entries ``velocity`` or ``temperature`` needs to be present.
   The program will skip initialization of a particular field if it cannot find
   an entry in the input file. This can be used to speed up the execution
   process if the user intends to initialize uniform velocity throughout the
   domain within Nalu.

``mesh_local_refinement``
-------------------------

This task creates an *error indicator field* that can be used to locally refine
the mesh using Percept. This is used to refine the wind farm simulation mesh
around the turbines to capture the wakes with the desired resolution while
performing the ABL simulations with a coarser mesh resolution.

.. code-block:: yaml

   mesh_local_refinement:
     fluid_parts: [fluid_part]
     write_percept_files: true
     percept_file_prefix: adapt
     search_tolerance: 11.0

     turbine_diameters: [ 15.0, 15.0 ]
     turbine_heights: [ 50.0, 50.0 ]
     turbine_locations:
       - [ 200.0, 200.0, 0.0 ]
       - [ 230.0, 300.0, 0.0 ]
     orientation:
       type: wind_direction
       wind_direction: 225.0
     refinement_levels:
       - [ 7.0, 12.0, 7.0 ]
       - [ 5.0, 10.0, 5.0 ]
       - [ 3.0, 6.0, 3.0]
       - [ 1.5, 3.0, 1.2]

.. confval:: turbine_diameters

   A list of turbine diameters for the turbines in the wind farm.

.. confval:: turbine_heights

   The list of tower heights for the turbines in the wind farm.

.. confval:: turbine_locations

   The ``(x, y, z)`` coordinates of the turbine base in the wind farm.

.. confval:: orientation

   The orientation of the refinement boxes. Currently there is only one option
   available indicated by ``type`` parameter: ``wind_direction``. For this
   option, it expects the ``wind_direction`` variable to contain the compass
   direction in degrees.

.. confval:: refinement_levels

   A list of 3 parameters for each nested refinement zone. The three parameters
   are the distance upstream, distance downstream and the lateral distance.
   These parameters are non-dimensional and are internally scaled by the turbine
   diameters by the utility. The nested boxes must be specified with the largest
   box first and the subsequent sizes in descending order.

.. confval:: search_tolerance

   The tolerance parameter added when searching for elements enclosed by the
   refinement box. A value slightly larger than the coarsest mesh size is
   recommended.

.. confval:: refine_field_name

   The name of the ``error_indicator_field`` used when creating STK fields.
   Default is ``turbine_refinement_field``.

.. confval:: write_percept_files

   Boolean flag indicating whether input files for use with Percept is written
   out by this utility as part of the run. Default: ``true``.

.. confval:: percept_file_prefix

   The prefix used for the Percept input file name. The default value is
   ``adapt``. With the default file name and three levels of refinement, it will
   create three input files: ``adapt1.yaml``, ``adapt2.yaml``, and
   ``adapt3.yaml``.

``init_channel_fields``
-----------------------

This task initializes the velocity fields for channel flow simulations
based on the parameters provided by the user and writes it out to the
:confval:`output_db`. It is safe to run ``init_channel_fields`` in
parallel. A sample invocation is shown below

.. code-block:: yaml
   :linenos:

   init_channel_fields:
     fluid_parts: [Unspecified-2-HEX]

     velocity:
       Re_tau : 550
       viscosity : 0.0000157

.. confval:: fluid_parts

   A list of element block names where the velocity fields are to be
   initialized.

.. confval:: velocity

   A YAML dictionary containing two values: the friction Reynolds
   number, ``Re_tau``, and the kinematic ``viscosity``
   (:math:`m^2/s`).

``generate_planes``
-------------------

Generates horizontal planes of nodesets at given heights that are used for
sampling velocity and temperature fields during an ABL simulation. The resulting
spatial average at given heights is used within Nalu to determine the driving
pressure gradient necessary to achieve the desired ABL profile during the
simulation. This task is capable of running in parallel.

The horizontal extent of the sampling plane can be either prescribed manually,
or the program will use the bounding box of the input mesh. Note that the latter
approach only works if the mesh boundaries are oriented along the major axes.
The extent and orientation of the sampling plane is controlled using the
``boundary_type`` option in the input file.

.. confval:: boundary_type

   Flag indicating how the program should estimate the horizontal extents of the
   sampling plane when generating nodesets. Currently, two options are supported:

   ==============  ===========================================================
   Type            Description
   ==============  ===========================================================
   bounding_box    Automatically estimate based on bounding box of the mesh
   quad_vertices   Use user-provided ``vertices``
   ==============  ===========================================================

   This flag is optional, and if it is not provided the program defaults to
   using the ``bounding_box`` approach to estimate horizontal extents.

.. confval:: fluid_part

   A list of element block names used to compute the extent using bounding box approach.

.. confval:: heights

   A list of vertical heights where the nodesets are generated.

.. confval:: part_name_format

   A ``printf`` style string that takes one floating point argument ``%f``
   representing the height of the plane. For example, if the user desires to
   generate nodesets at 70m and 90m respectively and desires to name the plane
   ``zh_070`` and ``zh_090`` respectively, this can be achieved by setting
   ``part_name_format: zh_%03.0f``.

.. confval:: dx, dy

   Uniform resolutions in the x- and y-directions when generating nodesets. Used
   only when :confval:`boundary_type` is set to ``bounding_box``.

.. confval:: nx, ny

   Number of subdivisions of along the two axes of the quadrilateral provided.
   Given 4 points, ``nx`` will divide segments ``1-2`` and ``3-4``, and ``ny``
   will divide segments ``2-3`` and ``4-1``. Used only when
   :confval:`boundary_type` is set to ``quad_vertices``.

.. confval:: vertices

   Used to provide the horizontal extents of the sampling plane to the utility. For example

   .. code-block:: yaml

      vertices:
        - [250.0, 0.0]      # Vertex 1 (S-W corner)
        - [500.0, -250.0]   # Vertex 2 (S-E corner)
        - [750.0, 0.0]      # Vertex 3 (N-E corner)
        - [500.0, 250.0]    # Vertex 4 (N-W corner)

Example using custom vertices
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml
   :linenos:

   generate_planes:
     boundary_type: quad_vertices       # Override default behavior
     fluid_part: Unspecified-2-hex      # Fluid part

     heights: [ 70.0 ]                  # Heights were sampling planes are generated
     part_name_format: "zplane_%06.1f"  # Name format for new nodesets
     nx: 25                             # X resolution
     ny: 25                             # Y resolution
     vertices:                          # Vertices of the quadrilateral
       - [250.0, 0.0]
       - [500.0, -250.0]
       - [750.0, 0.0]
       - [500.0, 250.0]


``create_bdy_io_mesh``
----------------------

Create an I/O transfer mesh containing the boundaries of a given ABL precursor
mesh. The I/O transfer mesh can be used with Nalu during the precursor runs to
dump inflow planes for use with a later wind farm LES simulation with
inflow/outflow boundaries. Unlike other utilities described in this section,
this utility creates a new mesh instead of adding to the database written out by
the :program:`nalu_preprocess` executable. It is safe to invoke this task in a
parallel MPI run.

.. confval:: output_db

   Name of the I/O transfer mesh where the boundary planes are written out. This
   argument is mandatory.

.. confval:: boundary_parts

   A list of boundary parts that are saved in the I/O mesh. The names in the
   list must correspond to the names of the sidesets in the given ABL mesh.


``rotate_mesh``
---------------

Rotates the mesh given angle, origin, and axis using quaternion rotations.

.. confval:: mesh_parts

   A list of element block names that must be rotated.

.. confval:: angle

   The rotation angle in degrees.

.. confval:: origin

   An (x, y, z) coordinate for mesh rotation.

.. confval:: axis

   A unit vector about which the mesh is rotated.

.. code-block:: yaml
   :linenos:

   rotate_mesh:
     mesh_parts:
       - unspecified-2-hex

     angle: 30.0
     origin: [500.0, 0.0, 0.0]
     axis: [0.0, 0.0, 1.0]
