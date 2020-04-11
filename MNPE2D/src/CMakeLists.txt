cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(MNPE2D)

enable_language(Fortran)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
    set(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -Wall")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-g -traceback -ffpe-trap=zero,invalid,overflow,underflow")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ip -xHOST")
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
    set(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -Wall")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -ffpe-trap=zero,invalid,overflow,underflow")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native")
endif()

add_executable(
    MNPE2D
    Pemp2dbb.f
    Envprop1.f
    phsprop.f
    pesrc.f
    Zsgen.f
    fft.f
    ssi.f
    )