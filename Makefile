#================================================================================
#
# COMPILER SPECIFIC SETTINGS
#
#================================================================================
# uncomment this section for Intel Fortran
FC                 := ifort

# # uncomment this section for gfortran
#FC                 := gfortran
# Also uncomment these next two lines and set ICCLIBPATH to correct
# path on your system if using gfortran but with a netCDF library that
# was compiled with Intel C This will get rid of undefined refs to
# intel_fast_memset, etc. You need to find the right path for this...
#ICCLIBPATH = <INTELPATH>
#LIB = -L$(ICCLIBPATH) -lirc

#================================================================================
#
# LIBRARY SETTINGS
#
#================================================================================

# Set to the path that contains the netecdf lib and include directories on your system
NETCDFPATH         :=    /nopt/nrel/apps/netcdf/dist/netcdf

#================================================================================
#
# COMPILERS AND FLAGS
#
#================================================================================
LD                 := $(FC)
FCFLAGS            := -O3 -g
LDFLAGS            := -lnetcdf
INCLUDES           := -I$(NETCDFPATH)/include
LIBS               := -L$(NETCDFPATH)/lib $(LIB)

#================================================================================
#
# BUILD
#
#================================================================================
BINNAME            := wrftonalu

build : $(BINNAME)

module_dm.o : module_dm.F90
	$(FC) -o $@ -c $(FCFLAGS) $(INCLUDES) $<

module_constants.o : module_constants.F90
	$(FC) -o $@ -c $(FCFLAGS) $(INCLUDES) $<

module_str2int.o : module_str2int.F90
	$(FC) -o $@ -c $(FCFLAGS) $(INCLUDES) $<

module_ncderrcheck.o : module_ncderrcheck.F90
	$(FC) -o $@ -c $(FCFLAGS) $(INCLUDES) $<

module_utmdeg_converter.o : module_utmdeg_converter.F90 module_str2int.o module_constants.o
	$(FC) -o $@ -c $(FCFLAGS) $(INCLUDES) $<

module_exodus.o : module_exodus.F90 module_dm.o module_ncderrcheck.o module_utmdeg_converter.o
	$(FC) -o $@ -c $(FCFLAGS) $(INCLUDES) $<

wrftonalu.o : wrftonalu.F90 module_dm.o module_constants.o module_ncderrcheck.o module_exodus.o
	$(FC) -o $@ -c $(FCFLAGS) $(INCLUDES) $<

$(BINNAME) : module_dm.o module_constants.o module_str2int.o module_ncderrcheck.o module_utmdeg_converter.o module_exodus.o wrftonalu.o 
	$(LD) -o $@ $^ $(LIBS) $(LDFLAGS) $(LIB)

doc : 
	@echo "Generating documentation with Doxygen."
	doxygen Doxyfile

clean :
	@echo "Clean of objects, modules, doc, and executable"
	@-rm -rf *.o *.mod doc $(BINNAME)

.PHONY : clean build doc
