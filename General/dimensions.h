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

#ifndef DIMENSIONS_H
#define DIMENSIONS_H

#include "../MPI/Decomposition/decomposition.h"
#include "macro.h"
#include "structs.h"

/*!
 * @brief This class is used to store information on dimensions of the data.
 */
class Dimensions {
    IndicesIJ elts_glob;            // Total number of elements in the domain in each direction
    IndicesIJ elts_loc;             // Number of local elements in each direction
    double dx;                      // Grid step size along x axis
    double dy;                      // Grid step size along y axis
    IndicesBegEnd internal_range_i; // beg/end indices of the internal elements in i-th direction
    IndicesBegEnd internal_range_j; // beg/end indices of the internal elements in j-th direction
    IndicesIJ elts_loc_with_halo;   // Total number of elements in all directions (including halo elements)
    Decomposition decomp;           // Stores information on the domain decomposition
    IndicesIJ beg_ind_glob;         // Global indices that determine the very first cell on the current sub-domain

public:
    /*!
     * @brief Default constructor.
     */
    Dimensions() : elts_glob(1, 1),
                   elts_loc(1, 1),
                   internal_range_i(0, 0),
                   internal_range_j(0, 0),
                   elts_loc_with_halo(1, 1),
                   dx(1.), dy(1.)
                   { }

    /* ************************************************* */
    /* ***************** Main functions **************** */
    /* ************************************************* */
    /*!
     * @brief Set up the geometrical properties of the domain.
     * @param L - length of the domain.
     */
    inline void setupGeometry(double L) {
        dx = L / elts_glob.i;
        dy = L / elts_glob.j;
    }

    /*!
     * @brief Decompose the domain.
     * @param num_procs_i Number of processes in i-th direction.
     * @param num_procs_j Number of processes in j-th direction.
     * @return EXIT_FAILURE if decomposition fails, EXIT_SUCCESS otherwise.
     */
    inline int decompose(const IndicesIJ &num_procs) {
        int exit_code = decomp.decompose(num_procs, elts_glob, elts_loc, beg_ind_glob);
        if (exit_code == EXIT_SUCCESS)
            findInternalIndices();
        return exit_code;
    }

    /* ************************************************* */
    /* ******************** Getters ******************** */
    /* ************************************************* */
    /*!
     * \returns Range of internal indices in i-th direction.
     */
    inline IndicesBegEnd getInternalIndRangeI() const { return internal_range_i; }

    /*!
     * \returns Range of internal indices in j-th direction.
     */
    inline IndicesBegEnd getInternalIndRangeJ() const { return internal_range_j; }

    /*!
     * \returns Number of internal elements in each direction.
     */
    inline IndicesIJ getNumEltsLoc() const { return elts_loc; }

//    /*!
//     * \returns Number of internal elements in j-th direction.
//     */
//    inline int getNumEltsLocJ() const { return elts_loc.j; }

    /*!
     * \returns Total number of elements in each direction, including halo elements.
     */
    inline IndicesIJ getNumElts() const { return elts_loc_with_halo; }

//    /*!
//     * \returns Total number of elements in j-th direction, including halo elements.
//     */
//    inline int getNumEltsJ() const { return elts_loc_with_halo.j; }

    /*!
     * \returns Constant reference to the object of decomposition.
     */
    inline const Decomposition &getDecomposition() const { return decomp; }

    /*!
     * \returns Grid step size along the x-axis.
     */
    inline double getDx() const { return dx; }

    /*!
     * \returns Grid step size along the y-axis.
     */
    inline double getDy() const { return dy; }

    /*!
     * \returns Total number of elements in un-decomposed domain in each direction.
     */
    inline IndicesIJ getNumEltsGlob() const { return elts_glob; }

    inline IndicesIJ getBegIndicesGlob() const { return beg_ind_glob; }

    /* ************************************************* */
    /* ******************** Setters ******************** */
    /* ************************************************* */
    /*!
     * @brief Set total number of elements in un-decomposed domain in each direction.
     * @param _imax Global number of elements in j-th direction.
     */
    inline void setNumEltsGlob(const IndicesIJ num_elts) { elts_glob = num_elts; }

    /*!
     * @brief Set local number of elements in the domain in each direction.
     * @param _imax_loc Local number of elements in i-th direction, excluding halo elements.
     */
    inline void setNumEltsLoc(const IndicesIJ &num_elts) { elts_loc = num_elts; }

private:
    /*!
     * @brief Find the range of the internal indices, including halo cells.
     */
    void findInternalIndices();
};

#endif
