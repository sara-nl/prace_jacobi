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
 * @file helpers.cpp
 * @brief Contains definitions of miscellaneous functions.
 */

#include <sys/time.h>
#include <iostream>
#include <cstdlib>
#include <string>
#include "helpers.h"
#include "structs.h"
#include "../MPI/common.h"

using namespace std;

void Helpers::setDimensionsAndDecompose(int argc, char** argv, Dimensions &dims) {

    IndicesIJ elts_glob;    // Number of global cells in each direction
    IndicesIJ num_procs;    // Number of processes in each direction

    parseInput(argc, argv, elts_glob, num_procs);

    /* Decompose the domain and assign local Dimensions */
    dims.setNumEltsGlob(elts_glob);
    if (dims.decompose(num_procs) == EXIT_FAILURE) {
        terminateExecution();
    }

}

void Helpers::parseInput(int argc, char** argv, IndicesIJ &elts_glob, IndicesIJ &num_procs) {

    /* Assign the default values first. */
    elts_glob.i = elts_glob.j = 10;
    num_procs.i = num_procs.j = 1;

    if (argc > 1 && argc != 7) {
        terminateDueToParserFailure();
    }
    else if (argc == 7) {
        int position_s = 1;
        int position_p = 4;

        /* Check the first key and assume the position of the second one */
        if (string(argv[1]) == "-s") {
            position_s = 1;
            position_p = 4;
        }
        else if (string(argv[1]) == "-d") {
            position_s = 4;
            position_p = 1;
        }

        /* Check that all keys are correct */
        if (std::string(argv[position_s]) != "-s") {
            terminateDueToParserFailure();
        }
        if (std::string(argv[position_p]) != "-d") {
            terminateDueToParserFailure();
        }

        /* Read the values */
        elts_glob.i = atoi(argv[position_s + 1]);
        elts_glob.j = atoi(argv[position_s + 2]);

        num_procs.i = atoi(argv[position_p + 1]);
        num_procs.j = atoi(argv[position_p + 2]);
    }
}

void Helpers::terminateDueToParserFailure() {
    printByRoot("\nError! Incorrect arguments were passed to the command line.\n"
                "Use the following keys:\n"
                "  -s - set number of the grid cells in each direction (i j)\n"
                "  -d - set decomposition for each direction (i j)\n"
                "Example:\n"
                "  ./a.out -s 10 10 -d 1 1");
    terminateExecution();
}

double Helpers::tic() {

    double current_time = 0.0;
#ifdef USE_MPI
    /* Always synchronize the processes before starting the timer. */
    MPI_Barrier(MPI_COMM_WORLD);
    current_time = MPI_Wtime();
#else
    struct timeval tv;
    gettimeofday(&tv, NULL);
    current_time = (double)tv.tv_sec + 1.0e-6 * (double) tv.tv_usec;
#endif

    return current_time;
}

double Helpers::toc() {

    double current_time = 0.0;
#ifdef USE_MPI
    /* We do not need to sync the processes when the timer is stopped because we are
     * interested in which process reached this point last. */
    current_time = MPI_Wtime();
    MPI_Barrier(MPI_COMM_WORLD);
#else
    struct timeval tv;
    gettimeofday(&tv, NULL);
    current_time = (double)tv.tv_sec + 1.0e-6 * (double) tv.tv_usec;
#endif

    return current_time;
}