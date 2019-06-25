# NHGOMO
Non-hydrostatic Generalized Operator Modelling of the Ocean (NHGOMO) is a three-dimensional ocean model based on OpenArray (Available at https://github.com/hxmhuang/OpenArray_CXX) which is a simple operator library for the decoupling of ocean modelling and parallel computing. For more details about OpenArray, please see the paper (https://www.geosci-model-dev-discuss.net/gmd-2019-28/) and the simple user mannual of OpenArray (Located at https://github.com/hxmhuang/OpenArray_CXX/tree/master/doc)

NHGOMO is a non-hydrostatic extention for GOMO (https://github.com/hxmhuang/GOMO). The fundamental equations and algorithms of GOMO are derived from POM2k (Blumberg and Mellor, 1987). GOMO features bottom-following, free-surface, staggered Arakawa C-grid. To effectively evolve the rapid surface fluctuations, GOMO uses the mode-splitting algorithm to address the fast propagating surface gravity waves and slow propagating internal waves in barotropic (external) and baroclinic (internal) modes, respectively.

NHGOMO is composed of Fortran files (.F90), header files (.h), a single namelist file (.txt), and a makefile. You can run the ideal test--lock-exchange problem, the default input file lock_exchange_flow.nc is located at the directory ./bin/data. 


# Compile NHGOMO
To compile NHGOMO, the following software are required:
  1) OpenArray (You can follow the installation guide in the user mannual step by step).
  2) Fortran 90 or Fortran 95 compiler.
  3) gcc/g++ compiler version 6.1.0 or higher. 
  4) Intel icc/icpc compiler version 2017 or higher. 
  5) GNU make version 3.81 or higher. 
  6) Message Passing Interface (MPI) library. 
  7) Parallel NetCDF library version 1.7.0 or higher. 
  8) Armadillo, a C++ library for linear algebra & scientific computing, version 8.200.2 or higher. 
  9) Boost C++ Libraries, version 1.65.1 or higher.

After the installation of the required software is done, You can follow the basic steps to compile NHGOMO:

  1) Download NHGOMO from GitHub:
        git clone https://github.com/hxmhuang/NHGOMO.git;
  2) Set environment variables: edit ./env_set to specify path to the required software and library, then type: . env_set;
  3) Edit the makefile: set the path to the libopenarray.a and openarray.mod of OpenArray;
  4) Make: make main 

# Run NHGOMO
After compiling, the executable file ./bin/nh_gomo will be generated. Within the directory ./bin where nh_gomo and namelist file (config.txt) exist, type:
  ./nh_gomo 
or
  mpirun -np N ./nh_gomo

