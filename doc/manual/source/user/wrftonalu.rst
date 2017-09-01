.. _util_wrftonalu_exe:

``wrftonalu`` -- WRF to Nalu Convertor
======================================

This program converts WRF data to the `Nalu <https://github.com/NaluCFD/Nalu>`_
(Exodus II) data format. Exodus II is part of
`SEACAS <https://gsjaardema.github.io/seacas>`_ and one can find other utilities to work
with Exodus II files there. The objective is to provide Nalu with input WRF data
as boundary conditions.

This program was started as ``WRFTOOF``, a WRF to OpenFoam converter,
which was written by J. Michalakes and M. Churchfield. It was adapted
for converting to Nalu data by M. T. Henry de Frahan.

.. note::

   This utility is not built by default. The user must set
   :confval:`ENABLE_WRFTONALU` to ``ON`` during the :ref:`CMake configure phase
   <user_installation>`.

Command line invocation
-----------------------

.. code-block:: bash

   bash$ wrftonalu [options] wrfout

where :file:`wrfout` is the WRF data file used to generate inflow conditions for
the Nalu simulations. The user must provide the relevant boundary files in the
run directory named :file:`west.g`, :file:`east.g`, :file:`south.g`,
:file:`north.g`, :file:`lower.g`, and :file:`upper.g`. Only the boundaries where
inflow data is requires needs to exist. The interpolated WRF data is written out
to files with extension ``*.nc`` for the corresponding grid files for use with
Nalu. The following optional parameters can be supplied to customize the
behavior of :program:`wrftonalu`.

.. program:: wrftonalu

.. option:: -startdate

   Date string of the form ``YYYY-mm-dd_hh_mm_ss`` or ``YYYY-mm-dd_hh:mm:ss``

.. option:: -offset

   Number of seconds to start Exodus directory naming (default: 0)

.. option:: -coord_offset lat lon

   Latitude and longitude of origin for Exodus mesh. Default: center of WRF data.

.. option:: -ic

   Populate initial conditions as well as boundary conditions.

.. option:: -qwall

   Generate temperature flux for the terrain (lower) BC file.
