#================================================================================
#
# Configuration
#
#================================================================================
SRCDIR             := src
BINDIR             := bin
TARGET             := wrftonalu

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
FCFLAGS            := -g -C
FCFLAGS            := -O3 -g
LDFLAGS            := -lnetcdf
INCLUDES           := -I$(NETCDFPATH)/include
LIBS               := -L$(NETCDFPATH)/lib $(LIB)

# Make sure the mod files go in the oject directory
ifeq ($(FC),ifort)
  MODULE_OUTPUT_FLAG := -module
else ifeq ($(FC),gfortran)
  MODULE_OUTPUT_FLAG := -J
endif

#================================================================================
#
# BUILD
#
#================================================================================

# Path where objects and dependencies will be put
ODIR?=objects

# files to compile
srcf90:=$(subst ./,,$(basename $(wildcard $(SRCDIR)/*.F90)))
objs:=$(foreach file,$(srcf90),$(ODIR)/$(file).o)

all: $(BINDIR)/$(TARGET)

$(ODIR)/%.o : %.F90
	@mkdir -p $(@D)
	@echo "Compilation of $< "
	$(FC) -o $@ -c $(FCFLAGS) $(INCLUDES) $< $(MODULE_OUTPUT_FLAG) $(@D)

$(BINDIR)/$(TARGET) : $(objs)
	@echo "Link of $@"
	@mkdir -p $(@D)
	$(LD) -o $@ $(FCFLAGS) $(objs) $(LIBS) $(LDFLAGS) 
	@echo "Successful compilation"

doc : 
	@echo "Generating documentation with Doxygen."
	doxygen Doxyfile

clean :
	@echo "Clean of $(ODIR), doc, and $(BINDIR)"
	@-rm -rf $(ODIR) $(BINDIR) doc

cleanall :
	@-rm -rv $(ODIR) $(BINDIR) doc

.PHONY : clean cleanall all doc
