.. _util_nalu_preprocess_exe:


``nalu_preprocess`` -- Nalu Preprocessing Utilities
===================================================

This utility loads an input mesh and performs various pre-processing *tasks* so
that the resulting output database can be used in a wind LES simulation.
Currently, the following *tasks* have been implemented within this utility.

====================  ===========================================================
Task type             Description
====================  ===========================================================
``init_abl_fields``   Initialize ABL velocity and temperature fields
``generate_planes``   Generate horizontal sampling planes for ``dp/dx`` forcing
``rotate_mesh``       Rotate mesh
====================  ===========================================================

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

   Filename where the pre-processed results database is output, e.g., ``ablNeutralPrecursor.g``

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



``init_abl_fields``
-------------------

This task initializes the vertical velocity and temperature profiles for use
with an ABL precursor simulations based on the parameters provided by the user
and writes it out to the :confval:`output_db`. It is safe to run
``init_abl_fields`` in parallel. A sample invocation is shown below

.. literalinclude:: files/nalu_preprocess.yaml
   :language: yaml
   :linenos:
   :lines: 27-44

.. confval:: fluid_parts

   A list of element block names where the velocity and/or temperature fields
   are to be initialized.

.. confval:: temperature

   A YAML dictionary containing two arrays: ``heights`` and the corresponding
   ``values`` at those heights. The data must be provided in SI units. No
   conversion is performed within the code.

.. confval:: velocity

   A YAML dictionary containing two arrays: ``heights`` and the corresponding
   ``values`` at those heights. The data must be provided in SI units. No
   conversion is performed within the code. The values in this case are two
   dimensional lists of shape ``[nheights, 3]`` where ``nheights`` is the length
   of the `heights` array provided.

.. note::

   Only one of the entries ``velocity`` or ``temperature`` needs to be present.
   The program will skip initialization of a particular field if it cannot find
   an entry in the input file. This can be used to speed up the execution
   process if the user intends to initialize uniform velocity throughout the
   domain within Nalu.

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


``calc_ndtw2d``
---------------

.. deprecated:: 0.1.0

   The implementation uses a brute-force method to compute the nearest wall
   distance and as such is unsitable for production use. Use only on small
   two-dimensional meshes.

Calculate the nearest distance to wall (NDTW) for 2-D airfoil meshes.


.. code-block:: yaml
   :linenos:

   calc_ndtw2d:
     fluid_parts:
       - Unspecified-2-QUAD
       - Unspecified-3-QUAD

     wall_parts:
       - airfoil
