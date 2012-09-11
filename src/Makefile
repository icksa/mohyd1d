MOAB_DIR=$(HOME)/Downloads/moab/moab-opt-parallel
MOAB_DIR_DEBUG=$(HOME)/Downloads/moab/moab-debug-parallel
HDF5_DIR=/usr/local/hdf5-1.8.8-parallel
NCDF_DIR=/usr/local/netcdf-4.1.3-serial

CXXFLAGS=-I$(MOAB_DIR)/include \
         -I$(HDF5_DIR)/include \
         -I$(NCDF_DIR)/include \
         -O2 -pedantic -g -Wall

CXXFLAGS_PROF=-pg -g

CXXFLAGS_DEBUG=-I$(MOAB_DIR_DEBUG)/include \
               -I$(HDF5_DIR)/include \
               -I$(NCDF_DIR)/include \
               -g -pedantic -Wall

LDFLAGS=-L$(MOAB_DIR)/lib -lMOAB \
        -L$(HDF5_DIR)/lib -lhdf5 \
        $(NCDF_DIR)/lib/libnetcdf.a -lz \
        -O2 -g

LDFLAGS_DEBUG=-L$(MOAB_DIR)/lib -lMOAB \
        -L$(HDF5_DIR)/lib -lhdf5 \
        $(NCDF_DIR)/lib/libnetcdf.a -lz \
	-g

LDFLAGS_PROF=-pg -g

all:
	mpicxx $(CXXFLAGS) -c main.cxx
	mpicxx main.o $(LDFLAGS) -o mohyd1d

debug:
	mpicxx $(CXXFLAGS_DEBUG) -c main.cxx
	mpicxx main.o $(LDFLAGS_DEBUG)  -o mohyd1d_debug

prof:
	mpicxx $(CXXFLAGS) $(CXXFLAGS_PROF) -c main.cxx
	mpicxx main.o $(LDFLAGS) $(LDFLAGS_PROF) -o mohyd1d_prof

asm:
	mpicxx $(CXXFLAGS) -S main.cxx

clean:
	rm -rf main.o
	rm -rf mohyd1d 
	rm -rf mohyd1d_debug
	rm -rf mohyd1d_prof
	rm -rf *.vtk
