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
 * @file vector.h
 * @brief Contains declaration of the \e Vector class.
 */

#ifndef VECTOR_H
#define VECTOR_H

#ifdef USE_MPI
#include <mpi.h>
#endif
#include <iostream>
#include <vector>
#include "../General/dimensions.h"
#include "../General/structs.h"
#include "matrix.h"

using namespace std;

/*!
 * @brief Structure of the on-border elements IDs.
 */
struct vector_ngb_ids {
    std::vector<int> east;
    std::vector<int> west;
    std::vector<int> south;
    std::vector<int> north;
};

/*!
 * @class Vector
 * @brief Represents dense vector.
 * The vector can be local or distributed.
 */
class Vector : public Matrix {
    // already have # of real elements and # of halo elements
    Neighbors halo_chunk_size;              // Number of halo elements in each
                                            // direction
    Neighbors halo_chunk_start_index;       // Starting index of halo elements
                                            // in each direction
    vector_ngb_ids on_boarder_ids;          // IndicesBegEnd of on-boarder elements
                                            // that should be sent to neighboring
                                            // processes

public:
    /*!
     * @brief Default constructor
     */
    Vector() {
        rows = 0;
        cols = 1;
        _loc_elts = 0;
        _halo_elts = 0;
     }

    /*!
     * @brief Allocate memory for the vector.
     * @note The memory is allocated but not initialized.
     * This method also identifies on-border elements which values should be
     * communicated to the neighboring processes.
     * @param in_dims [in] Dimensions of the numerical problem.
     */
    void resize(Dimensions const &in_dims);

    /*!
     * @brief Return a reference to the element at specified index (row)
     * @param row [in] Row.
     */
    inline double &operator()(int row) {
        return data[row];
    }

    /*!
     * @brief Transfer the data from the real cells of the local process to
     *        the halo cells of the remote process.
     */
    void exchangeRealHalo();

private:
    /*!
     * @brief Calculate chunk size and its starting index for the halo elements.
     * @param num_elts [in] Number of halo elements
     * @param _halo_start_index [in/out] Counter
     * @param _chunk_size [out] Chunk size
     * @param _chunk_start_index [out] Index of the first element in chunk
     */
    void associateChunkData(const int num_elts, int &_halo_start_index,
                            int &_chunk_size, int &_chunk_start_index);


#ifdef USE_MPI
    void sendAndReceive(vector<double> &snd_buf, vector<double> &rcv_buf, 
                        int pid, int tag, MPI_Request *request);

    void packSndBuffer(vector<double> &snd_buf, std::vector<int> &on_boarder_ids, 
                    int halo_chunk_size);

    void unpackRcvBuffer(vector<double> &rcv_buf, int halo_chunk_start_index, 
                        int halo_chunk_size, int ngb_pid);

    void communicateToNgb(vector<double> &snd_buf, vector<double> &rcv_buf,
                        vector<int> &on_boarder_ids, int halo_chunk_size,
                        int ngb_pid, int tag, MPI_Request *request);

    void waitForAll(Neighbors &ngb_pid,
                    MPI_Request *request_w, MPI_Request *request_e,
                    MPI_Request *request_s, MPI_Request *request_n);
#endif
};

#endif
