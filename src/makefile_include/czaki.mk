# Czaki laptop configuration (gnu-x86)
SYSTYPE = gnu-x86

ZOLTANLIB=-L/home/czaki/libs/zoltan/3.82/lib -lzoltan
SPRNGLIB=-L/usr/lib -lsprng
GMPLIB=-lgmp
HYPRELIB=-L/home/czaki/libs/hypre/2.10/lib -lHYPRE
ZOLTANINC=-I/home/czaki/libs/zoltan/3.82/include
SPRNGINC=-I/usr/include/sprng
HYPREINC= -I/home/czaki/libs/hypre/2.10/include
MPIINC=-I/usr/include/mpi
