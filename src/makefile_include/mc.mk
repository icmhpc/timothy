# MC laptop configuration
SYSTYPE = gnu-x86
ZOLTANLIB=-L/opt/zoltan/3.81/lib -lzoltan
SPRNGLIB=-L/usr/lib -lsprng
GMPLIB=-lgmp
HYPRELIB=-L/opt/hypre/2.9.0b/lib -lHYPRE -lblas -llapack
ZOLTANINC=-I/opt/zoltan/3.81/include
SPRNGINC=-I/usr/include/sprng
HYPREINC=-I/opt/hypre/2.9.0b/include
MPIINC=-I/usr/include/mpi
