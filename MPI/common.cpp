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
 * @file common.cpp
 * @brief Contains C-based definitions of the main MPI functions
 */

#include <iostream>
#include "common.h"
#include "../General/macro.h"

void findGlobalMin(double &value) {

#ifdef USE_MPI
    /*
     * Hint:
     *  - use MPI_IN_PLACE to "replace" the value
     */
    NOT_IMPLEMENTED
#endif
}

void findGlobalMax(double &value) {

#ifdef USE_MPI
    /*
     * Hint:
     *  - use MPI_IN_PLACE to "replace" the value
     */
    NOT_IMPLEMENTED
#endif
}

int getMyRank() {

    int my_rank = 0;
#ifdef USE_MPI
    NOT_IMPLEMENTED
#endif
    return my_rank;
}

int getNumProcs() {

    int num_procs = 1;
#ifdef USE_MPI
    NOT_IMPLEMENTED
#endif
    return num_procs;
}

void findGlobalSum(double &value) {
#ifdef USE_MPI
    /*
     * Hint:
     *  - use MPI_IN_PLACE to "replace" the value
     */
    NOT_IMPLEMENTED
#endif
}

void findGlobalSum(int &value) {
#ifdef USE_MPI
    /*
     * Hint:
     *  - use MPI_IN_PLACE to "replace" the value
     */
    NOT_IMPLEMENTED
#endif
}

void initialize(int argc, char** argv) {
#ifdef USE_MPI
    MPI_Init(&argc, &argv);
#endif
}

void finalize() {
#ifdef USE_MPI
    MPI_Finalize();
#endif
}

void printByRoot(const std::string& str) {
    int my_rank = getMyRank();
    if (my_rank == 0) {
        std::cout << str << "\n";
    }
}
