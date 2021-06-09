/*
 * Copyright (c) 2021 Maksim Masterov, SURF
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

/*!
 * @file io.cpp
 * @brief Contains definitions of methods from the \e IO class.
 */

#include <fstream>
#include "io.h"

// Writes data into the file
void IO::writeFile(std::string file_name, Dimensions &dims, Field &T) {

    printByRoot("Writing results to file: " + file_name);

#ifndef USE_MPI
    ofstream out;

    out.open(file_name);

    for(int i = 0; i < T.numRows(); ++i) {
        for(int j = 0; j < T.numCols(); ++j) {
            out << dims.getDx() * i + 0.5 * dims.getDx() << " "
                << dims.getDy() * j + 0.5 * dims.getDy() << " "
                << T(i, j) << "\n";
        }
    }

    out.close();
#else
    // MPI_File mpi_file;
    // int out_case = IO_BY_ROOT;

    // /*
    //  * Be sure you rewrite the file!
    //  * Just delete the old one
    //  */
    // NOT_IMPLEMENTED

    // switch (out_case) {
    //     case IO_BY_ROOT: default:
    //         // Using MPI_Gatherv() and writing by the root process
    //         writeByRoot(mpi_file, dims, T);
    //         break;

    //     case IO_BY_COLLECTIVE:
    //         // Performing parallel IO
    //         writeByAll(mpi_file, dims, T);
    //         break;
    // }

    // MPI_File_close(&mpi_file);
#endif
}

#ifdef USE_MPI
void IO::generateGrid(Dimensions &dims, Field &T, vector<double> &grid_1D) {

    // Assemble 1D array from the grid.
    int start_i = dims.getBegIndicesGlob().i;
    int start_j = dims.getBegIndicesGlob().j;
    int counter = 0;

    for(int i = 0; i < T.numRows(); ++i) {
        for(int j = 0; j < T.numCols(); ++j) {
            grid_1D[counter++] = dims.getDx() * (i + start_i) + 0.5 * dims.getDx();
            grid_1D[counter++] = dims.getDy() * (j + start_j) + 0.5 * dims.getDy();
        }
    }
}

void IO::convertTo1D(Field &field_2D, vector<double> &field_1D) {

    int counter = 0;
    for(int i = 0; i < field_2D.numRows(); ++i) {
        for(int j = 0; j < field_2D.numCols(); ++j) {
            field_1D[counter++] = field_2D(i, j);
        }
    }
}

void IO::writeByRoot(MPI_File &mpi_file, Dimensions &dims, Field &T) {

    MPI_Status status;
    int num_glob_elts = T.getDimensions().getNumEltsGlob().i    // global number of elements in the field
                        * T.getDimensions().getNumEltsGlob().j;
    vector<double> buffer_T_rcv(num_glob_elts, 0.0);            // buffer for receiving the filed
    vector<double> buffer_grid_rcv(2 * num_glob_elts, 0.0);     // buffer for receiving the grid
    vector<double> buffer_wrt(3 * num_glob_elts, 0.0);          // writing buffer
    vector<double> T_1D(T.size(), 0.0);                         // field stored as a 1D array
    vector<double> grid_1D(2 * T.size(), 0.0);                  // grid stored as a 1D array (x0,y0;x1,y1;...)
    vector<int> all_sizes(getNumProcs(), 0);                    // vector of all sizes from all processes
    vector<int> displacement(getNumProcs(), 0);                 // displacement used to indicate the "partitioning" of data
    int root_pid = 0;                                           // PID of the process that will perform writing
    int local_size = 0;                                         // number of local elements

    // Assemble 1D array from the field.
    convertTo1D(T, T_1D);

    // Assemble 1D array from the grid.
    generateGrid(dims, T, grid_1D);

    // For the field:
    /* ******************************************************** */
    // Gather the number of elements from each process
    // According to the standard: "The root process receives the messages and stores them in rank order."
    local_size = T.size();
    NOT_IMPLEMENTED

    // Calculate the displacement
    for(int i = 0; i < all_sizes.size(); ++i) {
        for(int j = 0; j < i; ++j) {
            displacement[i] += all_sizes[j];
        }
    }

    // Gather the data from the distributed field
    NOT_IMPLEMENTED
    /* ******************************************************** */

    // Do the same for the grid. Note, that we have two coordinates
    /* ******************************************************** */
    local_size = 2 * T.size();
    NOT_IMPLEMENTED

    // Calculate the displacement
    for(int i = 0; i < displacement.size(); ++i) {
        displacement[i] = 0;
    }
    for(int i = 0; i < all_sizes.size(); ++i) {
        for(int j = 0; j < i; ++j) {
            displacement[i] += all_sizes[j];
        }
    }

    // Gather the data from the distributed field
    NOT_IMPLEMENTED
    /* ******************************************************** */

    // Write the gathered data by the root process
    if (getMyRank() == root_pid) {
        int counter_grid = 0;
        int counter = 0;
        // Assemble the data in a single buffer
        for(int i = 0; i < T.getDimensions().getNumEltsGlob().i; ++i) {
            for(int j = 0; j < T.getDimensions().getNumEltsGlob().j; ++j) {
                buffer_wrt[counter + 0] = buffer_grid_rcv[counter_grid + 0];
                buffer_wrt[counter + 1] = buffer_grid_rcv[counter_grid + 1];
                buffer_wrt[counter + 2] = buffer_T_rcv[j + T.getDimensions().getNumEltsGlob().j * i];
                counter_grid += 2;
                counter += 3;
            }
        }
        // Write to the file
        NOT_IMPLEMENTED
    }
}

void IO::writeByAll(MPI_File &mpi_file, Dimensions &dims, Field &T) {

    /* Prepare new data structure */
    struct CombinedType {
        double grid[2];
        double temp;
    };
    vector<CombinedType> combined_data(T.size());

    /* Combine grid points and temperature values in a single Array of Structures (AoS) */
    int start_i = dims.getBegIndicesGlob().i;
    int start_j = dims.getBegIndicesGlob().j;
    for(int i = 0; i < T.numRows(); ++i) {
        for(int j = 0; j < T.numCols(); ++j) {
            combined_data[j + i * T.numCols()].grid[0] = dims.getDx() * (i + start_i) + 0.5 * dims.getDx();
            combined_data[j + i * T.numCols()].grid[1] = dims.getDy() * (j + start_j) + 0.5 * dims.getDy();
            combined_data[j + i * T.numCols()].temp = T(i, j);
        }
    }

    /* Prepare new MPI type */
    MPI_Datatype combined_data_type;
    MPI_Aint base_address;
    MPI_Aint displacements[2];
    MPI_Datatype types[2] = {MPI_DOUBLE, MPI_DOUBLE};
    int lengths[2] = {2, 1};
    struct CombinedType dummy_struct;

    /* Find the addresses and calculate the displacement */
    MPI_Get_address(&dummy_struct, &base_address);
    MPI_Get_address(&dummy_struct.grid[0], &displacements[0]);
    MPI_Get_address(&dummy_struct.temp, &displacements[1]);
    displacements[0] = MPI_Aint_diff(displacements[0], base_address);
    displacements[1] = MPI_Aint_diff(displacements[1], base_address);

    /* Create a new type of structure and commit it */
    NOT_IMPLEMENTED

    /* Prepare subarray type */
    int glob_sizes[2] = {T.getDimensions().getNumEltsGlob().i, T.getDimensions().getNumEltsGlob().j};
    int loc_sizes[2] = {T.getDimensions().getNumEltsLoc().i, T.getDimensions().getNumEltsLoc().j};
    int start_inds[2] = {T.getDimensions().getBegIndicesGlob().i, T.getDimensions().getBegIndicesGlob().j};
    MPI_Datatype subarray_type;

    /* Based on the created type of structure, create a subarray and commit it */
    NOT_IMPLEMENTED

    /* Set vile view and write data */
    MPI_Status status;
    /* Set the file view using the created types of structure and subarray */
    NOT_IMPLEMENTED

    /* Write to the file using "regular" and collective write */
    NOT_IMPLEMENTED
}
#endif
