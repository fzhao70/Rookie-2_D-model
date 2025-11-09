# Makefile for Rookie 2D Model
# Author: Fanghe Zhao
# A simple 2D atmospheric model for teaching

# ============================================================================
# Compiler settings
# ============================================================================
# Default compiler (can be overridden: make FC=ifort)
FC = gfortran

# Compiler flags
FFLAGS = -O2 -fdefault-real-8 -fopenmp
DEBUGFLAGS = -g -Wall -Wextra -fbounds-check -fbacktrace

# NetCDF settings - adjust paths if needed
NETCDF_INC = $(shell nf-config --fflags 2>/dev/null || echo "-I/usr/include")
NETCDF_LIB = $(shell nf-config --flibs 2>/dev/null || echo "-lnetcdff -lnetcdf")

# Combine flags
FCFLAGS = $(FFLAGS) $(NETCDF_INC)
LDFLAGS = $(NETCDF_LIB)

# ============================================================================
# Target executables
# ============================================================================
TARGETS = simple_model geo_wind

# ============================================================================
# Source files
# ============================================================================
SIMPLE_MODEL_SRC = simple_model_type.f90
GEO_WIND_SRC = geo_wind.f90

# Object files
SIMPLE_MODEL_OBJ = $(SIMPLE_MODEL_SRC:.f90=.o)
GEO_WIND_OBJ = $(GEO_WIND_SRC:.f90=.o)

# ============================================================================
# Rules
# ============================================================================

.PHONY: all clean debug help

# Default target
all: $(TARGETS)

# Simple model executable
simple_model: $(SIMPLE_MODEL_OBJ)
	@echo "Linking simple_model..."
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)
	@echo "Successfully built simple_model"

# Geostrophic wind executable
geo_wind: $(GEO_WIND_OBJ)
	@echo "Linking geo_wind..."
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)
	@echo "Successfully built geo_wind"

# Pattern rule for object files
%.o: %.f90
	@echo "Compiling $<..."
	$(FC) $(FCFLAGS) -c $< -o $@

# Debug build
debug: FFLAGS = $(DEBUGFLAGS)
debug: clean all
	@echo "Debug build complete"

# Clean build artifacts
clean:
	@echo "Cleaning build artifacts..."
	rm -f $(TARGETS) *.o *.mod *.MOD
	@echo "Clean complete"

# Clean all including output files
cleanall: clean
	@echo "Cleaning output files..."
	rm -f *.nc
	@echo "All clean"

# Help target
help:
	@echo "Rookie 2D Model - Makefile Help"
	@echo "================================"
	@echo ""
	@echo "Available targets:"
	@echo "  all         - Build all executables (default)"
	@echo "  simple_model - Build the simple 2D model"
	@echo "  geo_wind    - Build the geostrophic wind program"
	@echo "  debug       - Build with debug flags"
	@echo "  clean       - Remove executables and object files"
	@echo "  cleanall    - Remove all generated files including NetCDF output"
	@echo "  help        - Show this help message"
	@echo ""
	@echo "Usage examples:"
	@echo "  make              # Build all programs"
	@echo "  make simple_model # Build only simple_model"
	@echo "  make FC=ifort     # Build with Intel Fortran compiler"
	@echo "  make debug        # Build with debug flags"
	@echo ""
	@echo "Requirements:"
	@echo "  - Fortran compiler (gfortran or ifort)"
	@echo "  - NetCDF Fortran library"
	@echo ""

# ============================================================================
# Dependencies
# ============================================================================
# The simple_model_type.f90 contains multiple modules, so we need to ensure
# they are compiled in the correct order. Since they're all in one file,
# a single compilation handles this.

$(SIMPLE_MODEL_OBJ): $(SIMPLE_MODEL_SRC)
$(GEO_WIND_OBJ): $(GEO_WIND_SRC)
