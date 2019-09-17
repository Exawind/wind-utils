.. _util_slice_mesh_exe:

``slice_mesh`` -- Sampling plane generation
===========================================

The ``slice_mesh`` executable can be used to generate sampling planes that can
be used with I/O transfer interface of Nalu-Wind for extract subsets of data
from wind farm simulations.

Command line invocation
-----------------------

.. code-block:: bash

   bash$ slice_mesh -i slice_mesh.yaml

   Slice Mesh Generation Utility
   Input file: slice_mesh.yaml
   Loading slice inputs...
   Initializing slices...
   Slice: Registering parts to meta data:
     -  turbine1_1
     -  turbine1_2
   Slice: Registering parts to meta data:
     -  turbine2_1
     -  turbine2_2
   Generating slices for: turbine1
   Creating nodes... 10% 20% 30% 40% 50% 60% 70% 80% 90% 100%
   Creating elements... 10% 20% 30% 40% 50% 60% 70% 80% 90% 100%
   Generating coordinate field
    - turbine1_1
    - turbine1_2
   Generating slices for: turbine2
   Creating nodes... 10% 20% 30% 40% 50% 60% 70% 80% 90% 100%
   Creating elements... 10% 20% 30% 40% 50% 60% 70% 80% 90% 100%
   Generating coordinate field
    - turbine2_1
    - turbine2_2
   Writing mesh to file: sampling_planes.exo

   Memory usage: Avg:   10.957 MB; Min:   10.957 MB; Max:   10.957 MB

.. program:: slice_mesh

.. option:: -i, --input-file

   YAML input file to be processed. Default: ``slice_mesh.yaml``
