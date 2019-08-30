.. _util_boxturb_exe:

``boxturb`` -- Turbulence box utility
=====================================

The ``boxturb`` executable is used to convert binary turbulence files into
NetCDF format that can be read during Nalu-Wind simulations. In addition to
conversion, it allows the user to apply divergence correction and scaling the
different components through the input file.

Command line invocation
-----------------------

.. code-block:: bash

   bash$ boxturb -i boxturb.yaml

   Nalu Turbulent File Processing Utility
   Input file: boxturb.yaml
   Begin loading WindSim turbulence data
   	Loading file: sim1u.bin
   	Loading file: sim1v.bin
   	Loading file: sim1w.bin
   Begin output in NetCDF format: turbulence.nc
   NetCDF file written successfully: turbulence.nc

.. program:: boxturb

.. option:: -i, --input-file

   YAML inout file that contains inputs for the executable. Default: `boxturb.yaml`

Sample input file
-----------------

.. code-block:: yaml
   :linenos:

   boxturb:
     data_format: windsim
     output: turbulence.nc

     box_dims: [1024, 128, 128]
     box_len: [2400.0, 160.0, 160.0]

     bin_filenames:
       - sim1u.bin
       - sim1v.bin
       - sim1w.bin

     correct_divergence: yes

     solver_settings:
       method: pfmg
       preconditioner: none
       max_iterations: 200
       tolerance: 1.0e-8
       print_level: 1
       log_level: 1

     # Scaling factor
     apply_scaling: yes
     scale_type: default
     scaling_factors: [1.0, 0.7, 0.3]
