#!/bin/bash

compiler='g++'
extra_flags=(-std=c++1y)

if [ $# -eq 0 ]
then
    echo "The compilation type in not specified. Please, use 'omp', 'mpi', 'gpu' or 'hybrid' "
    echo "to compile with OpenMP, MPI or hybrid parallelism, respectively."
    exit 1
elif [ $1 = "mpi" ]
then
    echo "Compiling with 'mpicxx'..."
    compiler='mpicxx'
    extra_flags=(-DUSE_MPI -lmpi)
elif [ $1 = "omp" ]
then
    echo "Compiling with 'g++'..."
    extra_flags=(-fopenmp -fopenmp-simd -fopt-info-vec-optimized -march=native)
elif [ $1 = "hybrid" ]
then
    echo "Compiling in a hybrid mode with 'mpicxx'..."
    compiler='mpicxx'
    extra_flags=(-fopenmp -DUSE_MPI -lmpi)
elif [ $1 = "gpu" ]
then
    echo "Compiling with 'g++' and support for the OpenMP offloading..."
    extra_flags=(-fopenmp -foffload=nvptx-none='-misa=sm_35 -Ofast -lm' -DUSE_GPU)
else
    echo "Incorrect compilation type is specified. Please, use 'omp', 'mpi' or 'hybrid' "
    echo "to compile with OpenMP, MPI or hybrid parallelism, respectively."
    exit 1
fi

if [ $# -eq 2 ]
then
    if [ $2 = "test" ]
    then
        echo "Compiling for testing..."
        extra_flags+='-DTEST '
    else
        echo "Incorrect second argument."
    fi
fi

$compiler \
    "${extra_flags[@]}" \
    -g3 -O2 --std=c++11 \
    IO/io.cpp \
    General/helpers.cpp \
    Solver/solver.cpp \
    Solver/solver_gpu.cpp \
    System/system.cpp \
    General/dimensions.cpp \
    main.cpp \
    MPI/common.cpp \
    MPI/Decomposition/decomposition.cpp \
    DataTypes/matrix.cpp \
    DataTypes/vector.cpp \
    DataTypes/field.cpp \
    Tests/utests.cpp
