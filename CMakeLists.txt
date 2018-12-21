CMAKE_MINIMUM_REQUIRED(VERSION 3.0)
SET(CMAKE_Fortran_COMPILER ifort)

PROJECT(mkl_perfm_test)

ENABLE_LANGUAGE(Fortran)

SET(CMAKE_Fortran_FLAGS "-r8 -g -align -mkl=parallel -mcmodel=medium -heap-arrays 0 -traceback")
SET(CMAKE_Fortran_FLAGS_RELEASE "-O3 -qopenmp -qopenmp-simd -unroll -xHost \
    -qopt-report2 -qopt-report-phase=vec -qopt-report-phase=loop")
SET(CMAKE_Fortran_FLAGS_DEBUG "-O0 -W1 -CB")

ADD_SUBDIRECTORY(mkl_libs)
set(utils_src utils.f90)

#axpy
set(axpy_src BLAS1/axpy.f90 BLAS1/m_axpy.f90)
ADD_EXECUTABLE(axpy ${axpy_src} ${utils_src})
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/mkl_libs)
TARGET_LINK_LIBRARIES(axpy mkl_spblas)

#gather
set(gather_src Gather/gather.f90 Gather/m_gather.f90)
add_executable(gather ${gather_src} ${utils_src})
include_directories(${CMAKE_BINARY_DIR}/mkl_libs)
target_link_libraries(gather mkl_spblas)