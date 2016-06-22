#* **************************************************************************
# * This file is part of Timothy
# *
# * Copyright (c) 2014/15 Maciej Cytowski
# * Copyright (c) 2014/15 ICM, University of Warsaw, Poland
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# *
# * *************************************************************************/


ifeq ($(wildcard makefile_include/default.mk),)
    $(error Please create file deafult.mk in src/makefile_include/ directory)
endif 

include makefile_include/default.mk

all: $(EXEC) $(VALIDATOR)
	



#############################################################################
# Prepare configuration for your system 
#
#ifeq ($(SYSTYPE),HOSTNAME_OF_YOUR_SYSTEM)
#CC = mpicc
#OPT = -O2 -fopenmp  -ffast-math -ftree-vectorize -finline-functions -Winline -g
#DEFS =
#INCLUDES = -I/usr/include/mpi -I/usr/include/sprng -I. -I/opt/zoltan/3.81/include -I/opt/sprng/include -I/opt/hypre/2.9.0b/include
#LIBS = -lm -L/opt/zoltan/3.81/lib -lzoltan -L/opt/sprng/lib -lsprng -lgmp -L/opt/hypre/2.9.0b/lib -lHYPRE -lblas -llapack
#endif
#############################################################################

ifeq ($(SYSTYPE),okeanos)
CC = cc
OPT = -O2 -O2 -qopenmp
DEFS = -DPROF
INCLUDES = -I. $(ZOLTANINC) $(SPRNGINC) $(HYPREINC) $(MPIINC)
LIBS = -lm $(ZOLTANLIB) $(SPRNGLIB) $(GMPLIB) $(HYPRELIB)
endif


# IBM BG/Q configuration
ifeq ($(SYSTYPE),ibm-bgq)
CC = mpixlc_r
OPT = -O3 -qstrict -qtune=qp -qarch=qp -qsimd=auto -qsmp=omp -g -qinline -qkeyword=inline
DEFS = -DPROF
INCLUDES = -I. $(ZOLTANINC) $(SPRNGINC) $(HYPREINC) $(MPIINC)
LIBS = -lm $(ZOLTANLIB) $(SPRNGLIB) $(GMPLIB) $(HYPRELIB) -L/opt/ibmmath/essl/5.1/lib64 -lesslbg -L/opt/ibmcmp/xlf/bg/14.1/lib64 -lxlf90_r -lxlsmp
endif

# General GNU x86 configuration
ifeq ($(SYSTYPE),gnu-x86)
CC = mpicc
OPT = -O2 -fopenmp  -ffast-math -ftree-vectorize -finline-functions -Winline -Wall -Wextra -g # -Wconversion
DEFS = #-DDEBUG
INCLUDES = -I. $(ZOLTANINC) $(SPRNGINC) $(HYPREINC) $(MPIINC)
LIBS = -lm $(ZOLTANLIB) $(SPRNGLIB) $(GMPLIB) $(HYPRELIB)
endif

# Intel compiler, Intel MIC architecture
ifeq ($(SYSTYPE),intel-mic)
CC = mpiicc 
OPT = -O3 -openmp 
DEFS = #-DDEBUG
INCLUDES = -I. $(ZOLTANINC) $(SPRNGINC) $(HYPREINC) $(MPIINC)
LIBS = -lm $(ZOLTANLIB) $(SPRNGLIB) $(GMPLIB) $(HYPRELIB)
endif

# IBM Power775 configuration
ifeq ($(SYSTYPE),ibm-p775)
CC = mpcc_r -q64
OPT = -O3 -qinline -qunroll -qhot -qsimd=auto -qtune=pwr7 -qarch=pwr7 -qsmp=omp -qcpluscmt -DPROF -qkeyword=inline #-DDEBUG
INCLUDES = -I. $(ZOLTANINC) $(SPRNGINC) $(HYPREINC) $(MPIINC)
LIBS = -lmassv -lmass -lm $(ZOLTANLIB) $(SPRNGLIB) $(GMPLIB) $(HYPRELIB) -lessl -lxlf90_r -lxlsmp
endif

INCLUDES = -I. $(ZOLTANINC) $(SPRNGINC) $(HYPREINC) $(MPIINC)


EXEC = ../timothy
VALIDATOR = ../validator



OPTIONS = $(OPT) $(DEFS)

OBJS = main.o fake_io.o utils.o cells.o domdec.o octree.o potential.o comm.o stats.o random.o init.o grid.o interp.o fields.o chemf.o  compute.o

INCL =  cells.h  chemf.h  comm.h  compute.h  domdec.h  fields.h  global.h  grid.h  init.h  inline.h  interp.h  io.h  octree.h  potential.h  random.h  stats.h   utils.h

CFLAGS = $(OPTIONS) $(INCLUDES) 

$(EXEC): $(OBJS)
	$(CC) $(OPTIONS) $(OBJS) $(LIBS) -o $(EXEC)

$(OBJS): $(INCL)

VALIDATOR_FLAGS = -O0 -Winline -Wall -Wextra -g
VALIDATOR_CXX_FLAGS = -std=c++0x -std=c++11
VALIDATOR_OBJS = validator/validator.o ini_manipulator/ini_manipulator.o validator.o

VALIDATOR_INLC = global.h ini_parser/ini_parser.h

#CXX = clang++

$(VALIDATOR): $(VALIDATOR_OBJS)
	$(CXX) $(VALIDATOR_FLAGS) $(VALIDATOR_CXX_FLAGS) $(VALIDATOR_OBJS) -o $(VALIDATOR)

validator.o: validator.cpp ini_manipulator/ini_manipulator.h
	$(CXX) -c $(VALIDATOR_FLAGS) $(VALIDATOR_CXX_FLAGS) $< -o $@

validator/validator.o: validator/validator.cpp ini_manipulator/ini_manipulator.h
	$(CXX) -c $(VALIDATOR_FLAGS) $(VALIDATOR_CXX_FLAGS) $< -o $@

ini_manipulator/ini_manipulator.o: ini_manipulator/ini_manipulator.cpp ini_manipulator/ini_manipulator.h
	$(CXX) -c $(VALIDATOR_FLAGS) $(VALIDATOR_CXX_FLAGS) $< -o $@


.PHONY: clean
clean:
	rm -f $(OBJS); rm -f $(EXEC); rm -f $(VALIDATOR_OBJS); rm -f $(VALIDATOR)
clean_validator:
	rm -f $(VALIDATOR_OBJS); rm -f $(VALIDATOR)


doc:
	doxygen Doxyfile


