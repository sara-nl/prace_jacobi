/*
 * Copyright (c) 2024 Maksim Masterov, SURF
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
 * @file matrix.cpp
 * @brief Contains definition of methods from the \e Matrix class.
 */

#include "matrix.h"

void Matrix::resize(Dimensions const &in_dims) {

    int imax_loc = in_dims.getNumEltsLoc().i;
    int jmax_loc = in_dims.getNumEltsLoc().j;

    dims = in_dims;

    _loc_elts = dims.getNumEltsLoc().i * dims.getNumEltsLoc().j;
    _halo_elts = countHaloElts(dims);;

    rows = _loc_elts;
    cols = _loc_elts + _halo_elts;

    data.resize(rows * cols);
}

int Matrix::countHaloElts(Dimensions const &in_dims) {
    
    int imax_loc = dims.getNumEltsLoc().i;
    int jmax_loc = dims.getNumEltsLoc().j;
    int halo_elts = 0;

    if (in_dims.getDecomposition().getNgbPid().east != EMPTY)
        halo_elts += jmax_loc;

    if (in_dims.getDecomposition().getNgbPid().west != EMPTY)
        halo_elts += jmax_loc;

    if (in_dims.getDecomposition().getNgbPid().south != EMPTY)
        halo_elts += imax_loc;

    if (in_dims.getDecomposition().getNgbPid().north != EMPTY)
        halo_elts += imax_loc;

    return halo_elts;
}

void Matrix::print() {

    int procs = getNumProcs();
    int my_rank = getMyRank();

#ifdef USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    for(int pid = 0; pid < procs; ++pid) {
        if (pid == my_rank) {
            cout << "pid: " << pid << endl;
            for(int i = 0; i < rows; ++i) {
                for(int j = 0; j < cols; ++j) {
                    cout << data[j + i * cols] << " ";
                }
                cout << "\n";
            }
            cout << endl;
        }
#ifdef USE_MPI
        MPI_Barrier(MPI_COMM_WORLD);
#endif
    }
}
