#!/bin/bash

# Copyright (c) 2024 Maksim Masterov, SURF
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

compiler="icpc"
extra_flags=(-std=c++1y)


# ####################################### #
# Check version of the loaded MPI library #
# ####################################### #
function _check_mpi_version() {
    local mpi_compiler=""
    local mpi_type=""
    if ( ! command -v mpirun &> /dev/null )
    then
        echo "Error: The MPI library is not loaded (Abort)..." >&2
        echo "       Use 'module load 2023' to load the software environment" >&2
        echo "       and 'module av mpi/' to see available MPI libraries." >&2
        exit 1
    else
        mpi_type="$(mpirun --version | grep Open)"
        if [[ ! -z $mpi_type ]]
        then
            mpi_type="OpenMPI"
            echo "Loaded MPI library: OpenMPI" >&2
        else
            mpi_type="$(mpirun --version | grep Intel)"
            if [[ ! -z $mpi_type ]]
            then
                mpi_type="IntelMPI"
                echo "Loaded MPI library: Intel MPI" >&2
            else
                mpi_type="$(mpirun --version | grep HYDRA)"
                if [[ ! -z $mpi_type ]]
                then
                    mpi_type="MPICH"
                    echo "Loaded MPI library: MPICH" >&2
                else
                    echo "Error: The MPI library was found, but the version is unknown (Abort)..." >&2
                    exit 1
                fi
            fi
        fi
    fi

    if [[ $mpi_type == "OpenMPI" ]] || [[ $mpi_type == "MPICH" ]]
    then 
        mpi_compiler="mpicxx"
    else 
        mpi_compiler="mpiicpc"
    fi

    echo "$mpi_compiler"
}


# ####################################### #
#     Check name of the C++ compiler      #
# ####################################### #
function _check_cpp_compiler() {
    local cpp_compiler=""
    if ( ! command -v icpc &> /dev/null )
    then
        echo "Using GNU compiler..." >&2
        cpp_compiler="g++"
    else
        echo "Using Intel compiler..." >&2
        cpp_compiler="icpc"
    fi

    echo "$cpp_compiler"
}


# ####################################### #
#        Analyze input parameters         #
# ####################################### #
if [ $# -eq 0 ]
then
    echo "The compilation type in not specified. Please, use 'omp', 'mpi', 'hybrid', or 'gpu' " >&2
    echo "to compile with OpenMP, MPI, hybrid, or GPU parallelism, respectively." >&2
    exit 1
elif [ $1 = "mpi" ]
then
    read compiler < <( _check_mpi_version ) || exit 1
    echo "Compiling with '$compiler'..." >&2
    extra_flags=(-DUSE_MPI -lmpi)
elif [ $1 = "omp" ]
then
    read compiler < <( _check_cpp_compiler ) || exit 1
    echo "Compiling with '$compiler'..." >&2
    if [[ $compiler ==  "icpc" ]]
    then
        extra_flags+=(-qopenmp)
    else
        extra_flags+=(-fopenmp)
    fi
elif [ $1 = "hybrid" ]
then
    read compiler < <( _check_mpi_version ) || exit 1
    echo "Compiling in a hybrid mode with '$compiler'..." >&2
    extra_flags=(-DUSE_MPI -lmpi)
    if [[ $compiler ==  "mpiicpc" ]]
    then
        extra_flags+=(-qopenmp)
    else
        extra_flags+=(-fopenmp)
    fi
elif [ $1 = "gpu" ]
then
    read compiler < <( _check_cpp_compiler ) || exit 1
    echo "Compiling with $compiler' and support for the OpenMP offloading..." >&2
    extra_flags=(-fopenmp -foffload=nvptx-none='-misa=sm_35 -Ofast -lm')
else
    echo "Error: Incorrect compilation type is specified. Please, use 'omp', 'mpi' or 'hybrid' " >&2
    echo "       to compile with OpenMP, MPI or hybrid parallelism, respectively." >&2
    exit 1
fi

if [ $# -eq 2 ]
then
    if [ $2 = "test" ]
    then
        echo "Compiling for testing..." >&2
        extra_flags+=(-DTEST)
    else
        echo "Error: Incorrect second argument." >&2
        exit 1
    fi
fi


# ####################################### #
#            Compile the code             #
# ####################################### #
$compiler \
    "${extra_flags[@]}" \
    -g -O2 \
    IO/io.cpp \
    General/helpers.cpp \
    Solver/solver.cpp \
    System/system.cpp \
    General/dimensions.cpp \
    main.cpp \
    MPI/common.cpp \
    MPI/Decomposition/decomposition.cpp \
    DataTypes/matrix.cpp \
    DataTypes/vector.cpp \
    DataTypes/field.cpp \
    Tests/utests.cpp
