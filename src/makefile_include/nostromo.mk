# Nostromo (ibm-bgq) system configuration
SYSTYPE = ibm-bgq
ZOLTANLIB=-L/opt/zoltan/3.8/lib -lzoltan
SPRNGLIB=-L/opt/sprng/2.0/lib -lsprng
GMPLIB=-L/opt/gmp/5.1.1/lib -lgmp
HYPRELIB=-L/opt/hypre/2.9.0b_bigint_omp/lib -lHYPRE
ZOLTANINC=-I/opt/zoltan/3.8/include
SPRNGINC=-I/opt/sprng/2.0/include
HYPREINC=-I/opt/hypre/2.9.0b_bigint_omp/include
MPIINC=
