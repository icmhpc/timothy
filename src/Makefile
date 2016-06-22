ifeq ($(wildcard makefile_include/default.mk),)
    $(error Please create file deafult.mk in src/makefile_include/ directory)
endif
include makefile_include/default.mk

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

# okeanos - ICM's Cray XC40 configuration
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
OPT = -O2 -fopenmp  -ffast-math -ftree-vectorize -finline-functions -Winline -Wall -Wextra -g
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

all: $(EXEC) $(VALIDATOR)

OPTIONS = $(OPT) $(DEFS)

OBJS = main.o io.o utils.o cells.o domdec.o octree.o potential.o comm.o stats.o random.o init.o grid.o interp.o fields.o tempf.o chemf.o vessel.o bone.o compute.o

INCL = bone.h  cells.h  chemf.h  comm.h  compute.h  domdec.h  fields.h  global.h  grid.h  init.h  inline.h  interp.h  io.h  octree.h  potential.h  random.h  stats.h  tempf.h  utils.h  vessel.h

CFLAGS = $(OPTIONS) $(INCLUDES) 

$(EXEC): $(OBJS)
	$(CC) $(OPTIONS) $(OBJS) $(LIBS) -o $(EXEC)

$(OBJS): $(INCL)

VALIDATOR_FLAGS = -O0 -Winline -Wall -Wextra -g
VALIDATOR_CXX_FLAGS = -std=c++0x -std=c++11
VALIDATOR_OBJS = validator.o ini_parser/ini_parser.o validator/validator.o

VALIDATOR_INLC = global.h ini_parser/ini_parser.h

#CXX = clang++

$(VALIDATOR): $(VALIDATOR_OBJS)
	$(CXX) $(VALIDATOR_FLAGS) $(VALIDATOR_CXX_FLAGS) $(VALIDATOR_OBJS) -o $(VALIDATOR)

validator.o: validator.cpp ini_parser/ini_parser.h
	$(CXX) -c $(VALIDATOR_FLAGS) $(VALIDATOR_CXX_FLAGS) $< -o $@

validator/validator.o: validator/validator.cpp ini_parser/ini_parser.h
	$(CXX) -c $(VALIDATOR_FLAGS) $(VALIDATOR_CXX_FLAGS) $< -o $@

ini_parser/ini_parser.o: ini_parser/ini_parser.c ini_parser/ini_parser.h
	cc -c $(VALIDATOR_FLAGS) $< -o $@

.PHONY: clean
clean:
	rm -f $(OBJS); rm -f $(EXEC); rm -f $(VALIDATOR_OBJS); rm -f $(VALIDATOR)
clean_validator:
	rm -f $(VALIDATOR_OBJS); rm -f $(VALIDATOR)


doc:
	doxygen Doxyfile
