CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(mkl_libs)

ENABLE_LANGUAGE(Fortran)

SET(CMAKE_Fortran_FLAGS "-r8 -g -align -mkl=parallel")
SET(CMAKE_Fortran_FLAGS_RELEASE "-O3 -xHost")
SET(CMAKE_Fortran_FLAGS_DEBUG "-O0 -W1 -CB")

ADD_LIBRARY(mkl_spblas STATIC ${MKL_ROOT}/include/mkl_spblas.f90)