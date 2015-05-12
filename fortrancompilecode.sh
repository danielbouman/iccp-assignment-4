#!/bin/bash
#Compile code:
gfortran -c timeevolution.f90
f2py3.4 -c --f90flags='-fPIC -O3 -funroll-loops -ffast-math -march=native' -m f90timeevolution timeevolution.f90 --fcompiler=gfortran
