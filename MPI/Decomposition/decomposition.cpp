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
 * @file decomposition.h
 * @brief Contains definition of methods from the class \e Decomposition.
 */

#include "decomposition.h"
#include <iostream>
#include <cmath>

using namespace std;

int Decomposition::getProcCoord(int &proc_ind_i, int &proc_ind_j) {

    int my_rank = getMyRank();

    if (num_subdomains.j == -1) {
        if (my_rank == 0) {
            cout << "Error! The function getProcCoord() can't be used, because the domain has not "
                    "beed decomposed yet." << endl;
        }
        return EXIT_FAILURE;
    }

    /*
     * Note, according to the standard, C and C++ always round down results of
     * the integer division!
     */
    proc_ind_i = my_rank / num_subdomains.j;                // index of the process
                                                            // in i-th direction
    proc_ind_j = my_rank - proc_ind_i * num_subdomains.j;   // index of the process
                                                            // in j-th direction

    return EXIT_SUCCESS;
}

int Decomposition::getProcInd(int proc_ind_i, int proc_ind_j) {

    return proc_ind_j + proc_ind_i * num_subdomains.j;
}

int Decomposition::decompose(const IndicesIJ num_procs, const IndicesIJ elts_glob,
                             IndicesIJ &elts_loc, IndicesIJ &beg_ind_glob) {

    int num_procs_avail = getNumProcs();

    /* Check if number of processes correspond to the decomposition size. */
    if (num_procs_avail != num_procs.i * num_procs.j) {
        printByRoot("The specified number of processes doesn't "
                    "match the available number of processes: "
                    + std::to_string(num_procs.i * num_procs.j)
                    + " vs. "
                    + std::to_string(num_procs_avail));
        /* Note, all processes should exit the function! */
        return EXIT_FAILURE;
    }

    /* Assign number of subdomains to local variables. */
    num_subdomains.i = num_procs.i;
    num_subdomains.j = num_procs.j;

    /* Assign local Dimensions */
    elts_loc.i = floor(elts_glob.i / num_subdomains.i);
    elts_loc.j = floor(elts_glob.j / num_subdomains.j);

    /*
     * Assume that all processes are enumerated in the "natural" order. for a 2d
     * decomposition among 9 processes the enumeration will look like:
     *   2 5 8
     *   1 4 7
     *   0 3 6
     * We need to ensure that processes at the most right column and the most top
     * row have correct number of elements. This is important in situations, when
     * the total number of elements in i-th or j-th directions are not divisible
     * by the number of specified processes in the same direction. For instance,
     * we have 3x3 decomposition (see above). Let's apply it to the domain of
     * 10x10 elements. Thus, processes 0,1,3,4 should have 3x3 elements each,
     * while processes 2 and 5 should have 3x4 elements, processes 6 and 7 should
     * have 4x3 elements, and process 8 should have 4x4 elements. Summing up, all
     * this processes will result in total number of 100 elements.
     */

    /* Get process "coordinates". Note: my_rank = proc_ind_j + proc_ind_i * nj. */
    int proc_ind_i = 0;     // index of the process in i-th direction.
    int proc_ind_j = 0;     // index of the process in j-th direction.
    if (getProcCoord(proc_ind_i, proc_ind_j) == EXIT_FAILURE)
        return EXIT_FAILURE;

    /*
     * Get the global indices that correspond to the very first (bottom-left) cell
     * of the current sub-domain.
     */
    beg_ind_glob.i = proc_ind_i * elts_loc.i;
    beg_ind_glob.j = proc_ind_j * elts_loc.j;

    /* Correct the number of cells */
    if ((elts_glob.i % num_subdomains.i) && (proc_ind_i == num_subdomains.i - 1)) {
        int remainder = elts_glob.i - elts_loc.i * (proc_ind_i + 1);
        elts_loc.i += remainder;
    }

    if ((elts_glob.j % num_subdomains.j) && (proc_ind_j == num_subdomains.j - 1)) {
        int remainder = elts_glob.j - elts_loc.j * (proc_ind_j + 1);
        elts_loc.j += remainder;
    }

    findNeighborsIds();
    checkForPhysicalBoundaries();

    return EXIT_SUCCESS;
}

int Decomposition::findNeighborsIds() {

    int my_rank = getMyRank();

    /* Get process "coordinates" */
    int proc_ind_i = 0;     // Index of the process in i-th direction.
    int proc_ind_j = 0;     // Index of the process in j-th direction.

    if (getProcCoord(proc_ind_i, proc_ind_j) == EXIT_FAILURE)
        return EXIT_FAILURE;

    ngb_pid.central = my_rank;

    /* Check the west neighbor. */
    if (proc_ind_i != 0)
        ngb_pid.west = getProcInd(proc_ind_i - 1, proc_ind_j);

    /* Check the east neighbor. */
    if (proc_ind_i != num_subdomains.i - 1)
        ngb_pid.east = getProcInd(proc_ind_i + 1, proc_ind_j);

    /* Check the south neighbor. */
    if (proc_ind_j != 0)
        ngb_pid.south = getProcInd(proc_ind_i, proc_ind_j - 1);

    /* Check the north neighbor. */
    if (proc_ind_j!= num_subdomains.j - 1)
        ngb_pid.north = getProcInd(proc_ind_i, proc_ind_j + 1);

    return EXIT_SUCCESS;
}

void Decomposition::checkForPhysicalBoundaries() {

    /*
     * Simply check if the neighboring subdomain exists. If not - there is a
     * physical boundary of that site.
     */
    if (ngb_pid.west == EMPTY)
        phys_bound.west = PHYS_BOUNDARY;
    if (ngb_pid.east == EMPTY)
        phys_bound.east = PHYS_BOUNDARY;
    if (ngb_pid.south == EMPTY)
        phys_bound.south = PHYS_BOUNDARY;
    if (ngb_pid.north == EMPTY)
        phys_bound.north = PHYS_BOUNDARY;
}
