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
 * @file main.cpp
 * @brief The starting point.
 * A simple MPI/OpenMP code that solves a 2D Poisson equation (no sources or
 * sinks) on a uniform structured grid using Jacobi method.
 */

#include "General/helpers.h"
#include "MPI/common.h"
#include "System/system.h"
#include "Solver/solver.h"
#include "IO/io.h"
#include "Tests/utests.h"

/*!
 * @brief Report elapsed time.
 * @param start Start time
 * @param end End time
 * @param message Message to be added to the report
 */
void reportElapsedTime(double start, double end, const string &message) {
    findGlobalMin(start);
    findGlobalMax(end);
    printByRoot("Elapsed time (" + message + "): " + std::to_string(end - start) + "s.");
}

/*!
 * @brief Run unit tests.
 */
int runTests() {

    Utests utests;
    return utests.runAll();
}

/*!
 * @brief Run a 2D Poisson problem.
 * @param argc [in] Number of command line arguments
 * @param argv [in] Vector of command line arguments
 */
void runProblem(int argc, char** argv) {

    Field T;                    // Temperature field
    Matrix A;                   // Matrix of the linear system
    Vector x, b;                // Vectors of unknowns and right hand side
    Dimensions dims;            // Dimensions of the problem
    Faces boundary_values;      // Boundary data
    System system;              // Object of the linear system
    Solver solver;              // Object of mathematical functions
    IO io;                      // Object for IO operations
    Helpers helpers;            // Object of auxiliary functions
    double elp_time[4] = {0};   // Elapsed time, [s]

    /* 
     * Check input from the command line and determine properties of the
     * numerical grid.
     */
    helpers.setDimensionsAndDecompose(argc, argv, dims);

    /* 
     * Set boundary values at walls. Note, all boundary conditions are
     * assumed to be of Dirichlet type.
     */
    boundary_values.east = 10.;
    boundary_values.west = 11.;
    boundary_values.south = 12.;
    boundary_values.north = 13.;

    /* Allocate memory for the distributed field, matrix and vectors. */
    system.allocateMemory(dims, T, A, x, b);

    /* Assemble the linear system. */
    system.assembleSystem(boundary_values, T, A, x, b);

    /* Solve the linear system. */
    elp_time[0] = helpers.tic();
    solver.solveJacobi(A, x, b);
    elp_time[1] = helpers.toc();

    /* Copy final solution back to the filed. */
    system.copySolution(x, T);

    /* Write results into the file. */
    elp_time[2] = helpers.tic();
    io.writeFile("output.dat", dims, T);
    elp_time[3] = helpers.toc();

    /* Report elapsed time. */
    reportElapsedTime(elp_time[0], elp_time[1], "Jacobi");
    reportElapsedTime(elp_time[2], elp_time[3], "IO");
}

int main(int argc, char** argv) {

    int exit_status = EXIT_SUCCESS;
    /*
     * Initialize the scope of MPI calls
     */
    initialize(argc, argv);

    /*
     * Choose to either run tests or to setup and run the problem
     */
#ifdef TEST
    exit_status = runTests();
#else
    runProblem(argc, argv);
#endif

    /* 
     * Finalize the scope of MPI calls
     */
    finalize();

    return exit_status;
}
