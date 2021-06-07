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
 * @brief Contains declaration of the class \e Decomposition.
 */

#ifndef DECOMPOSITION_H
#define DECOMPOSITION_H

#include "../common.h"
#include "../../General/macro.h"
#include "../../General/structs.h"
#include <vector>

/*!
 * @class Decomposition
 * @brief Responsible for the data decomposition in a 1D or 2D way.
 */
class Decomposition {
    IndicesIJ num_subdomains;   // Total number of subdomains in each direction
    Neighbors ngb_pid;          // Structure with indicators of the PIDs of the
                                // neighboring subdomains (EMPTY stands for "no neighbor")
    Neighbors phys_bound;       // Structure with indicators of the presence of
                                // the physical (real) boundary (PHYS_BOUNDARY
                                // stands for existing physical boundary)

public:
    /*!
     * @brief Default constructor.
     */
    Decomposition() : num_subdomains(1, 1) { }

    /*!
     * @brief Decompose the domain.
     * @param num_procs [in] Number of subdomains in each direction.
     * @param elts_glob [in] Global number of elements/cells in each direction.
     * @param elts_loc [out] Local number of elements/cells in each direction.
     * @param beg_ind_glob [out] Global indices of the very first cell.
     * @return Returns EXIT_SUCCESS on success and EXIT_FAILURE on error.
     */
    int decompose(const IndicesIJ num_procs, const IndicesIJ elts_glob,
                  IndicesIJ &elts_loc, IndicesIJ &beg_ind_glob);

    /*!
     * @brief Return a structure of the neighboring processes IDs.
     * If there is no neighboring process, the member of the structure is set
     * to EMPTY. Otherwise, the member is equal to the neighboring process ID.
     */
    inline const Neighbors &getNgbPid() const { return ngb_pid; }

    /*!
     * @brief Return a structure that indicates existence of the physical
     *        boundaries.
     * If member of the structure is equal to PHYS_BOUNDARY, the boundary exists.
     * Otherwise, the member is equal to EMPTY.
     */
    inline const Neighbors &getPhysBound() const { return phys_bound; }

private:
    /*!
     * @brief Evaluate process IDs of the neighboring sub-domains.
     * @return Returns EXIT_SUCCESS on success and EXIT_FAILURE on error.
     */
    int findNeighborsIds();

    /*!
     * @brief Determine if local sub-domain has any physical boundaries.
     */
    void checkForPhysicalBoundaries();

    /*!
     * @brief Calculate the coordinates of the sub-domain based on its rank.
     * @param proc_ind_i [out] Coordinate in i-th direction.
     * @param proc_ind_j [out] Coordinate in j-th direction.
     * @return Returns EXIT_SUCCESS on success and EXIT_FAILURE on error.
     */
    int getProcCoord(int &proc_ind_i, int &proc_ind_j);

    /*!
     * @brief Calculate the rank of the sub-domain based on its coordinates
     * @param proc_ind_i [in] Coordinate in i-th direction.
     * @param proc_ind_j [in] Coordinate in j-th direction.
     * @return Returns process ID.
     */
    int getProcInd(int proc_ind_i, int proc_ind_j);
};

#endif
