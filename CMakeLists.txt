CMAKE_MINIMUM_REQUIRED(VERSION 3.0)
SET(CMAKE_Fortran_COMPILER ifort)

PROJECT(mkl_perfm_test)

ENABLE_LANGUAGE(Fortran)

SET(SRCS main.f90 BLAS1/axpy.f90)
SET(CMAKE_Fortran_FLAGS "-r8 -g -align -mkl=parallel -heap-arrays 0 -traceback")
SET(CMAKE_Fortran_FLAGS_RELEASE "-O3 -qopenmp -unroll -xHost")
SET(CMAKE_Fortran_FLAGS_DEBUG "-O0 -W1 -CB")

ADD_SUBDIRECTORY(mkl_libs)

ADD_EXECUTABLE(${PROJECT_NAME} ${SRCS})
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/mkl_libs)
TARGET_LINK_LIBRARIES(${PROJECT_NAME} mkl_spblas)
