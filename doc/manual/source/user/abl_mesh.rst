.. _util_abl_mesh_exe:

``abl_mesh`` -- Block HEX Mesh Generation
=========================================

The ``abl_mesh`` executable can be used to generate structured mesh with HEX-8
elements in Exodus-II format. The interface is similar to OpenFOAM's
``blockMesh`` utility and can be used to generate simple meshes for ABL
simulations on flat terrain without resorting to commercial mesh generation
software, e.g., Pointwise.

Command line invocation
-----------------------

.. code-block:: bash

   bash$ abl_mesh -i abl_mesh.yaml

   Nalu ABL Mesh Generation Utility
   Input file: abl_mesh.yaml
   HexBlockMesh: Registering parts to meta data
   	Mesh block: fluid_part
   Num. nodes = 1331; Num elements = 1000
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
   Writing mesh to file: ablmesh.exo


.. program:: abl_mesh

.. option:: -i, --input-file

   YAML input file to be processed for mesh generation details. Default:
   ``nalu_abl_mesh.yaml``.

Input File Parameters
---------------------

The input file must contain a ``nalu_abl_mesh`` section that contains the input
parameters.A sample input file is shown below

.. code-block:: yaml
   :linenos:

   nalu_abl_mesh:
     output_db: ablmesh.exo

     spec_type: bounding_box

     vertices:
       - [0.0, 0.0, 0.0]
       - [10.0, 10.0, 10.0]

     mesh_dimensions: [10, 10, 10]


.. confval:: output_db [nalu_abl_mesh]

   The Exodus-II filename where the mesh is output. No default, must be provided
   by the user.

.. confval:: spec_type

   Specification type used to define the extents of the structured HEX mesh.
   This option is used to interpret the :confval:`vertices` read from the input
   file. Currently, two options are supported:

   =================  =======================================================
   Type               Description
   =================  =======================================================
   ``bounding_box``   Use axis aligned bounding box as domain boundaries
   ``vertices``       Use user provided vertices to define extents
   =================  =======================================================

.. confval:: vertices

   The coordinates specifying the extents of the computational domain. This
   entry is interpreted differently depending on the :confval:`spec_type`. If
   type is set to ``bounding_box`` then the code expects a list of two 3-D
   coordinate points describing bounding box to generate an axis aligned mesh.
   Otherwise, the code expects a list of 8 points describing the vertices of the
   trapezoidal prism.

.. confval:: mesh_dimensions

   Mesh resolution for the resulting structured HEX mesh along each direction.
   For a trapezoidal prism, the code will interpret the major axis along
   ``1-2``, ``1-4``, and ``1-5`` edges respectively.

.. confval:: fluid_part_name

   Name of the element block created with HEX-8 elements. Default value:
   ``fluid_part``.

.. confval:: ioss_8bit_ints

   Boolean flag that enables output of 8-bit ints when writing Exodus mesh.
   Default value: false.

Boundary names
~~~~~~~~~~~~~~

The user has the option to provide custom boundary names through the input file.
Use the boundary name input parameters to change the default parameters. If
these are not provided the default boundary names are described below:

======================  =====================
Boundary                Default sideset name
======================  =====================
``xmin_boundary_name``  ``west``
``xmax_boundary_name``  ``east``
``ymin_boundary_name``  ``south``
``ymax_boundary_name``  ``north``
``zmin_boundary_name``  ``terrain``
``zmax_boundary_name``  ``top``
======================  =====================

Mesh spacing
~~~~~~~~~~~~

Users can specify the mesh spacing to be applied in each direction by adding
additional sections (``x_spacing``, ``y_spacing``, and ``z_spacing``
respectively) to the input file. If no option is specified then a constant mesh
spacing is used in that direction.

========================== ===============================================
Available options          Implementation
========================== ===============================================
``constant_spacing``       :class:`~sierra::nalu::ConstantSpacing`
``geometric_stretching``   :class:`~sierra::nalu::GeometricStretching`
========================== ===============================================

**Example input file**

.. code-block:: yaml

   # Specifiy constant spacing in x direction (this is the default)
   x_spacing:
     spacing_type: constant_spacing

   # y direction has a mesh stretching factor
   y_spacing:
     spacing_type: geometric_stretching
     stretching_factor: 1.1

   # z direction has a mesh stretching factor in both directions
   z_spacing:
     spacing_type: geometric_stretching
     stretching_factor: 1.1
     bidirectional: true

Limitations
-----------

#. Does not support the ability to generate multiple blocks

#. Must be run on a single processor, running with multiple MPI ranks is currently
   unsupported.
