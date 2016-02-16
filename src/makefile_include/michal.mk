# Czaki laptop configuration (gnu-x86)
SYSTYPE = gnu-x86

ZOLTANLIB=-L/opt/zoltan/lib -lzoltan
ZOLTANINC=-I/opt/zoltan/include

SPRNGLIB=-L/usr/lib -lsprng
SPRNGINC=-I/usr/include/sprng

HYPRELIB=-L/opt/hypre/lib -lHYPRE
HYPREINC= -I/opt/hypre/include

MPIINC=-I/usr/include/mpi
GMPLIB=-lgmp
