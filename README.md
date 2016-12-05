# WRF to Nalu converter

This program converts WRF data to
the [Nalu](https://github.com/NaluCFD/Nalu) (Exodus II) data
format. Exodus II is part of [SEACAS](gsjaardema.github.io/seacas) and
one can find other utilities to work with Exodus II files there. The
objective is to provide Nalu with input WRF data as boundary
conditions.

This program was started as WRFTOOF, a WRF to OpenFoam converter,
which was written by J. Michalakes and M. Churchfield. It was adapted
for converting to Nalu data by M. T. Henry de Frahan.

## Building

The code  will compile and run with either gfortran or ifort and
needs to have the path to the NetCDF library and include files specified
in the [Makefile](Makefile). It's as easy as 

```{bash}
make
```

## Running

```{bash}
./wrftonalu wrfout
```
where `wrfout` is the WRF data file. This will generate Exodus
boundary condition files for Nalu.
