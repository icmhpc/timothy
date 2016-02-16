# Boreasz (ibm-p775) system configuration 
SYSTYPE = ibm-p775
ZOLTANLIB=-L/p7data/opt/zoltan/3.6/lib -lzoltan
SPRNGLIB=-L/p7data/opt/sprng/2.0/lib -lsprng
GMPLIB=-L/p7data/opt/gmp/5.0.5/lib -lgmp
HYPRELIB=-L/opt/hypre/2.9.0b/bigint_omp/lib -lHYPRE
ZOLTANINC=-I/p7data/opt/zoltan/3.6/include
SPRNGINC=-I/p7data/opt/sprng/2.0/include
HYPREINC=-I/opt/hypre/2.9.0b/bigint_omp/include
MPIINC=
