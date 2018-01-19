#!/bin/bash
# Shell script to compile CSimMDMV software V1.1
# This file was first created by Junwei Huang @ University of Toronto, Dec 2008.

mpiCC -c fftn.cpp
mpiCC -c -D_FILE_OFFSET_BITS=64 CSim3D3V.cpp
mpiCC -c -D_FILE_OFFSET_BITS=64 CSim3D2V.cpp
mpiCC -c -D_FILE_OFFSET_BITS=64 CSim3D1V.cpp

mpiCC -c -D_FILE_OFFSET_BITS=64 CSim2D3V.cpp
mpiCC -c -D_FILE_OFFSET_BITS=64 CSim2D2V.cpp
mpiCC -c -D_FILE_OFFSET_BITS=64 CSim2D1V.cpp

mpiCC -c -D_FILE_OFFSET_BITS=64 CSim3D3V_BENCH.cpp

mpiCC CSim3D3V.o fftn.o -o ../bin/CSim3D3V
mpiCC CSim3D2V.o fftn.o -o ../bin/CSim3D2V
mpiCC CSim3D1V.o fftn.o -o ../bin/CSim3D1V

mpiCC CSim2D3V.o fftn.o -o ../bin/CSim2D3V
mpiCC CSim2D2V.o fftn.o -o ../bin/CSim2D2V
mpiCC CSim2D1V.o fftn.o -o ../bin/CSim2D1V

mpiCC CSim3D3V_BENCH.o fftn.o -o ../bin/CSim3D3V_BENCH
rm *.o
