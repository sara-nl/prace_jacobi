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
 * @file system.h
 * @brief Contains definition of methods from the System class
 */

#include "../System/system.h"

void System::allocateMemory(Dimensions &dims, Field &T, Matrix &A,
                            Vector &x, Vector &b) {

    /* Allocate memory */
    x.resize(dims);
    b.resize(dims);
    A.resize(dims);
    T.resize(dims);

    /* Initialize data using first touch */
    for(int i = 0; i < A.numRows(); ++i) {
#pragma omp parallel for
        for(int j = 0; j < A.numCols(); ++j) {
            A(i, j) = 0.0;
        }
    }

#pragma omp parallel for
    for(int i = 0; i < x.numRows(); ++i) {
        x(i) = 0.0;
        b(i) = 0.0;
    }

    for(int i = 0; i < T.numRows(); ++i) {
#pragma omp parallel for
        for(int j = 0; j < T.numCols(); ++j) {
            T(i, j) = 0.0;
        }
    }
}

void System::assembleSystem(Faces &bondary_values, Field &T, Matrix &A,
                            Vector &x, Vector &b) {
    
    Faces coefficients;                                 // System coefficients
    int num_procs = 1;                                  // Total number of processes
    int my_rank = 0;                                    // Process rank
    IndicesBegEnd int_ind_i = T.getDimensions().getInternalIndRangeI(); // Pair of local begin/end
                                                        // IndicesBegEnd in i-th direction
    IndicesBegEnd int_ind_j = T.getDimensions().getInternalIndRangeJ(); // Pair of local begin/end
                                                        // IndicesBegEnd in j-th direction
    Dimensions dims = T.getDimensions();                // Problem Dimensions

    my_rank = getMyRank();
    num_procs = getNumProcs();

    /*
     * We are assembling a standard 5-point stencil using 2nd order central
     * difference scheme:
     *      [ 0 -1  0]
     *      [-1  4 -1]
     *      [ 0 -1  0]
     */
    coefficients.central = 4.;
    coefficients.east = -1.;
    coefficients.west = -1.;
    coefficients.south = -1.;
    coefficients.north = -1.;

    for(int i = int_ind_i.beg; i <= int_ind_i.end; ++i) {
        for(int j = int_ind_j.beg; j <= int_ind_j.end; ++j) {

            int row = T.getID(i, j);          // current row id
            int col = row;
            
            /* Central coefficient and corresponding LHS and RHS */
            A(row, row) = coefficients.central;
            b(row) = 0.0;
            x(row) = 0.0;
            
            /* Now, coefficients from the neighboring cells */
            /* On west */
            if (dims.getDecomposition().getPhysBound().west == PHYS_BOUNDARY && i == 0) {
                A(row, row) -= coefficients.west;
                b(row) -= 2. * coefficients.west * bondary_values.west;
            }
            else {
                col = T.getID(i - 1, j);
                A(row, col) = coefficients.west;
            }

            /* On east */
            if (dims.getDecomposition().getPhysBound().east == PHYS_BOUNDARY &&
                    i == T.getDimensions().getNumElts().i - 1) {
                A(row, row) -= coefficients.east;
                b(row) -= 2. * coefficients.east * bondary_values.east; 
            }
            else {
                col = T.getID(i + 1, j);
                A(row, col) = coefficients.east;
            }
            
            /* On south */
            if (dims.getDecomposition().getPhysBound().south == PHYS_BOUNDARY && j == 0) {
                A(row, row) -= coefficients.south;
                b(row) -= 2. * coefficients.south * bondary_values.south;
            }
            else {
                col = T.getID(i, j - 1);
                A(row, col) = coefficients.south;
            }

            /* On north */
            if (dims.getDecomposition().getPhysBound().north == PHYS_BOUNDARY &&
                    j == T.getDimensions().getNumElts().j - 1) {
                A(row, row) -= coefficients.north;
                b(row) -= 2. * coefficients.north * bondary_values.north;
            }
            else {
                col = T.getID(i, j + 1);
                A(row, col) = coefficients.north;
            }
        }
    }
}

void System::copySolution(Vector &x, Field &T) {

#pragma omp parallel for
    for(int i = 0; i < T.numRows(); ++i) {
        for(int j = 0; j < T.numCols(); ++j) {
            T(i, j) = x(j + i * T.numCols());
        }
    }
}
