#!/bin/bash

compiler="icpc"
extra_flags=(-std=c++11)

# ####################################### #
# Check version of the loaded MPI library #
# ####################################### #
function _check_mpi_version() {
    local mpi_compiler=""
    local mpi_type=""
    if ( ! command -v mpirun &> /dev/null )
    then
        echo "Error: The MPI library is not loaded (Abort)..." >&2
        echo "       Use 'module load 2021' to load the software environment" >&2
        echo "       and 'module av MPI' to see available MPI libraries." >&2
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
    echo "Compiling with 'g++' and support for the OpenMP offloading..." >&2
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
    -g -O3 \
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
