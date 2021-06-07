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

#include <cmath>
#include "utests.h"
#include "../MPI/common.h"
#include "../General/dimensions.h"
#include "../System/system.h"
#include "../Solver/solver.h"

void Utests::passed(const string name) {
    if (getMyRank() == 0)
        cout << name << " PASSED" << endl;
}

void Utests::failed(const string name) {
    if (getMyRank() == 0)
        cout << name << " FAILED" << endl;
}

int Utests::runAll() {
    int exit_status = 0;

    exit_status += decomposition1d();
    exit_status == EXIT_SUCCESS ? passed("1d decomposition                       ") :
                                  failed("1d decomposition                       ");

    exit_status += decomposition2d();
    exit_status == EXIT_SUCCESS ? passed("2d decomposition                       ") :
                                  failed("2d decomposition                       ");

    exit_status += matrixHalo1d();
    exit_status == EXIT_SUCCESS ? passed("matrix halo/real cells 1d decomposition") :
                                  failed("matrix halo/real cells 1d decomposition");

    exit_status += matrixHalo2d();
    exit_status == EXIT_SUCCESS ? passed("matrix halo/real cells 2d decomposition") :
                                  failed("matrix halo/real cells 2d decomposition");

    exit_status += vectorHalo1d();
    exit_status == EXIT_SUCCESS ? passed("vector halo/real cells 1d decomposition") :
                                  failed("vector halo/real cells 1d decomposition");

    exit_status += vectorHalo2d();
    exit_status == EXIT_SUCCESS ? passed("vector halo/real cells 2d decomposition") :
                                  failed("vector halo/real cells 2d decomposition");

    exit_status += fieldIDs2d();
    exit_status == EXIT_SUCCESS ? passed("enumeration of the field elements (2d) ") :
                                  failed("enumeration of the field elements (2d) ");

    exit_status += matrixAssembly2d();
    exit_status == EXIT_SUCCESS ? passed("matrix assembly (2d)                   ") :
                                  failed("matrix assembly (2d)                   ");

    exit_status += norm2d();
    exit_status == EXIT_SUCCESS ? passed("L2-norm (2d)                           ") :
                                  failed("L2-norm (2d)                           ");

    if (exit_status == 0)
        return EXIT_SUCCESS;
    else
        return EXIT_FAILURE;
}

int Utests::decomposition1d() {
    Dimensions dims;
    int check = EXIT_SUCCESS;
    int my_rank = getMyRank();
    IndicesIJ num_procs = {4, 1};

    dims.setNumEltsGlob({10, 10});

    dims.decompose(num_procs);

    if (dims.getNumEltsLoc().j != 10)
        check = EXIT_FAILURE;
    if (my_rank == 0 && dims.getNumEltsLoc().i != 2)
        check = EXIT_FAILURE;
    if (my_rank == 1 && dims.getNumEltsLoc().i != 2)
        check = EXIT_FAILURE;
    if (my_rank == 2 && dims.getNumEltsLoc().i != 2)
        check = EXIT_FAILURE;
    if (my_rank == 3 && dims.getNumEltsLoc().i != 4)
        check = EXIT_FAILURE;

    // This one is based on the assumtion that EXIT_SUCCESS is always 0
    findGlobalSum(check);
    return check > 0 ? EXIT_FAILURE : EXIT_SUCCESS;
}

int Utests::decomposition2d() {
    Dimensions dims;
    int check = EXIT_SUCCESS;
    int my_rank = getMyRank();
    IndicesIJ num_procs = {2, 2};

    dims.setNumEltsGlob({5, 5});

    dims.decompose(num_procs);

    if (my_rank == 0 && (dims.getNumEltsLoc().i != 2 || dims.getNumEltsLoc().j != 2))
        check = EXIT_FAILURE;
    if (my_rank == 1 && (dims.getNumEltsLoc().i != 2 || dims.getNumEltsLoc().j != 3))
        check = EXIT_FAILURE;
    if (my_rank == 2 && (dims.getNumEltsLoc().i != 3 || dims.getNumEltsLoc().j != 2))
        check = EXIT_FAILURE;
    if (my_rank == 3 && (dims.getNumEltsLoc().i != 3 || dims.getNumEltsLoc().j != 3))
        check = EXIT_FAILURE;

    // This one is based on the assumtion that EXIT_SUCCESS is always 0
    findGlobalSum(check);
    return check > 0 ? EXIT_FAILURE : EXIT_SUCCESS;
}

int Utests::matrixHalo1d() {
    Dimensions dims;
    int check = EXIT_SUCCESS;
    int my_rank = getMyRank();
    Matrix A;
    IndicesIJ num_procs = {4, 1};

    dims.setNumEltsGlob({5, 5});

    dims.decompose(num_procs);
    
    A.resize(dims);

    if (my_rank == 0 && (A.getLocElts() != 5 || A.getHaloElts() != 5 ||
                         A.numRows() != 5 || A.numCols() != 10))
        check = EXIT_FAILURE;
    if (my_rank == 1 && (A.getLocElts() != 5 || A.getHaloElts() != 10 ||
                         A.numRows() != 5 || A.numCols() != 15))
        check = EXIT_FAILURE;
    if (my_rank == 2 && (A.getLocElts() != 5 || A.getHaloElts() != 10 ||
                         A.numRows() != 5 || A.numCols() != 15))
        check = EXIT_FAILURE;
    if (my_rank == 3 && (A.getLocElts() != 10 || A.getHaloElts() != 5 ||
                         A.numRows() != 10 || A.numCols() != 15))
        check = EXIT_FAILURE;

    // This one is based on the assumtion that EXIT_SUCCESS is always 0
    findGlobalSum(check);
    return check > 0 ? EXIT_FAILURE : EXIT_SUCCESS;
}

int Utests::matrixHalo2d() {
    Dimensions dims;
    int check = EXIT_SUCCESS;
    int my_rank = getMyRank();
    Matrix A;
    IndicesIJ num_procs = {2, 2};

    dims.setNumEltsGlob({5, 5});

    dims.decompose(num_procs);
    
    A.resize(dims);

    if (my_rank == 0 && (A.getLocElts() != 4 || A.getHaloElts() != 4 || 
                         A.numRows() != 4 || A.numCols() != 8))
        check = EXIT_FAILURE;
    if (my_rank == 1 && (A.getLocElts() != 6 || A.getHaloElts() != 5 ||
                         A.numRows() != 6 || A.numCols() != 11))
        check = EXIT_FAILURE;
    if (my_rank == 2 && (A.getLocElts() != 6 || A.getHaloElts() != 5 ||
                         A.numRows() != 6 || A.numCols() != 11))
        check = EXIT_FAILURE;
    if (my_rank == 3 && (A.getLocElts() != 9 || A.getHaloElts() != 6 ||
                         A.numRows() != 9 || A.numCols() != 15))
        check = EXIT_FAILURE;

    // This one is based on the assumtion that EXIT_SUCCESS is always 0
    findGlobalSum(check);
    return check > 0 ? EXIT_FAILURE : EXIT_SUCCESS;
}

int Utests::vectorHalo1d() {
    Dimensions dims;
    int check = EXIT_SUCCESS;
    int my_rank = getMyRank();
    Vector x;
    IndicesIJ num_procs = {4, 1};

    dims.setNumEltsGlob({5, 5});

    dims.decompose(num_procs);
    
    x.resize(dims);

    if (my_rank == 0 && (x.getLocElts() != 5 || x.getHaloElts() != 5 ||
                         x.numRows() != 10 || x.numCols() != 1))
        check = EXIT_FAILURE;
    if (my_rank == 1 && (x.getLocElts() != 5 || x.getHaloElts() != 10 ||
                         x.numRows() != 15 || x.numCols() != 1))
        check = EXIT_FAILURE;
    if (my_rank == 2 && (x.getLocElts() != 5 || x.getHaloElts() != 10 ||
                         x.numRows() != 15 || x.numCols() != 1))
        check = EXIT_FAILURE;
    if (my_rank == 3 && (x.getLocElts() != 10 || x.getHaloElts() != 5 ||
                         x.numRows() != 15 || x.numCols() != 1))
        check = EXIT_FAILURE;

    // cout << my_rank << " " << x.getLocElts() << " " << x.getHaloElts() << " "
    //                 << x.numRows() << " " << x.numCols() << "\n";
    // This one is based on the assumtion that EXIT_SUCCESS is always 0
    findGlobalSum(check);
    return check > 0 ? EXIT_FAILURE : EXIT_SUCCESS;
}

int Utests::vectorHalo2d() {
    Dimensions dims;
    int check = EXIT_SUCCESS;
    int my_rank = getMyRank();
    Vector x;
    IndicesIJ num_procs = {2, 2};

    dims.setNumEltsGlob({5, 5});

    dims.decompose(num_procs);
    
    x.resize(dims);

    if (my_rank == 0 && (x.getLocElts() != 4 || x.getHaloElts() != 4 ||
                         x.numRows() != 8 || x.numCols() != 1))
        check = EXIT_FAILURE;
    if (my_rank == 1 && (x.getLocElts() != 6 || x.getHaloElts() != 5 ||
                         x.numRows() != 11 || x.numCols() != 1))
        check = EXIT_FAILURE;
    if (my_rank == 2 && (x.getLocElts() != 6 || x.getHaloElts() != 5 ||
                         x.numRows() != 11 || x.numCols() != 1))
        check = EXIT_FAILURE;
    if (my_rank == 3 && (x.getLocElts() != 9 || x.getHaloElts() != 6 ||
                         x.numRows() != 15 || x.numCols() != 1))
        check = EXIT_FAILURE;

    // This one is based on the assumtion that EXIT_SUCCESS is always 0
    findGlobalSum(check);
    return check > 0 ? EXIT_FAILURE : EXIT_SUCCESS;
}

int Utests::fieldIDs2d() {

    int ref_data[4][16] = {{0, 1, 4, 2, 3, 5, 6, 7, -1, -2, -2, -2, -2, -2, -2, -2},
                           {6, 0, 1, 2, 7, 3, 4, 5, -1, 8, 9, 10, -2, -2, -2, -2},
                           {6, 7, -1, 0, 1, 8, 2, 3, 9, 4, 5, 10, -2, -2, -2, -2},
                           {-1, 9, 10, 11, 12, 0, 1, 2, 13, 3, 4, 5, 14, 6, 7, 8}};
    int ref_data_size[4] = {9, 12, 12, 16};
    IndicesIJ num_procs = {2, 2};

    Dimensions dims;
    int check = EXIT_SUCCESS;
    int my_rank = getMyRank();
    Field field;

    dims.setNumEltsGlob({5, 5});

    dims.decompose(num_procs);
    
    field.resize(dims);

    if (ref_data_size[my_rank] != field.getIDs().size())
        check = EXIT_FAILURE;
    else {
        for(int i = 0; i < field.getIDs().size(); ++i) {
            if (ref_data[my_rank][i] != field.getIDs()[i])
                check = EXIT_FAILURE;
        }
    }

    // This one is based on the assumtion that EXIT_SUCCESS is always 0
    findGlobalSum(check);
    return check > 0 ? EXIT_FAILURE : EXIT_SUCCESS;
}

int Utests::matrixAssembly2d() {

    int ref_data[4][39] = {{6, -1, -1, -1,  5, -1, -1, -1,  5, -1, -1, -1, -1,  4, -1, -1, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2},
                           {5, -1, -1, -1, -1,  5, -1, -1, -1,  6, -1, -1,  4, -1, -1, -1, -1, -1,  4, -1, -1, -1, -1,  5, -1, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2},
                           {5, -1, -1, -1, -1,  4, -1, -1, -1, -1,  5, -1, -1, -1, -1,  4, -1, -1, -1,  6, -1, -1, -1,  5, -1, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2},
                           {4, -1, -1, -1, -1, -1,  4, -1, -1, -1, -1,  5, -1, -1, -1,  4, -1, -1, -1, -1, -1,  4, -1, -1, -1, -1,  5, -1, -1,  5, -1, -1, -1, -1,  5, -1, -1, -1,  6}};
    int ref_data_size[4] = {16, 25, 25, 39};
    IndicesIJ num_procs = {2, 2};

    System system;
    Dimensions dims;
    Faces boundary_values;
    int check = EXIT_SUCCESS;
    int my_rank = getMyRank();
    Field T;
    Matrix A;
    Vector x, b;
    int counter = 0;

    dims.setNumEltsGlob({5, 5});

    dims.decompose(num_procs);
    
    /* 
     * Setup boundary values at walls. Note, all boundary conditions are assumed to be 
     * of Dirichlet type
     */
    boundary_values.east = 10.;
    boundary_values.west = 11.;
    boundary_values.south = 12.;
    boundary_values.north = 13.;

    /* 
     * Assemble and solve the system
     */
    system.allocateMemory(dims, T, A, x, b);

    /* 
     * Assemble the linear system
     */
    system.assembleSystem(boundary_values, T, A, x, b);

    // check number of non-zero elements
    for(int i = 0; i < A.numRows(); ++i) {
        for(int j = 0; j < A.numCols(); ++j) {
            if (A(i, j) != 0) {
                ++counter;
            }
        }
    }
    
    if (ref_data_size[my_rank] != counter) {
        check = EXIT_FAILURE;
    }
    else {
        // check non-zero elements
        counter = 0;
        for(int i = 0; i < A.numRows(); ++i) {
            for(int j = 0; j < A.numCols(); ++j) {
                if (A(i, j) != 0) {
                    if (A(i, j) != ref_data[my_rank][counter])
                        check = EXIT_FAILURE;
                    ++counter;
                }
            }
        }
    }

    // This one is based on the assumtion that EXIT_SUCCESS is always 0
    findGlobalSum(check);
    return check > 0 ? EXIT_FAILURE : EXIT_SUCCESS;
}

int Utests::norm2d() {

    Solver solver;
    Dimensions dims;
    int check = EXIT_SUCCESS;
    int my_rank = getMyRank();
    Vector x;
    double answer = 36.3280883064331000;
    IndicesIJ num_procs = {2, 2};

    dims.setNumEltsGlob({5, 5});

    dims.decompose(num_procs);
    
    x.resize(dims);

    if (my_rank == 0) {
        x(0) = 3.1;
        x(1) = 4.8;
        x(2) = 9;
        x(3) = 3.5;
    }
    if (my_rank == 1) {
        x(0) = 1.1;
        x(1) = 7.4;
        x(2) = 3.3;
        x(3) = 2.9;
        x(4) = 5.5;
        x(5) = 11;
    }
    if (my_rank == 2) {
        x(0) = 9.2;
        x(1) = 4.4;
        x(2) = 1.4;
        x(3) = 3.9;
        x(4) = 7.3;
        x(5) = 8.4;
    }
    if (my_rank == 3) {
        x(0) = 8.6;
        x(1) = 1.1;
        x(2) = 9;
        x(3) = 6.5;
        x(4) = 9.9;
        x(5) = 16;
        x(6) = 7.7;
        x(7) = 8.9;
        x(8) = 5.6;
    }

    if (fabs(answer - solver.calculateNorm(x)) > 1e-12)
        check = EXIT_FAILURE;

    // This one is based on the assumtion that EXIT_SUCCESS is always 0
    findGlobalSum(check);
    return check > 0 ? EXIT_FAILURE : EXIT_SUCCESS;
}
