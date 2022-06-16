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
 * @file vector.cpp
 * @brief Contains definitions of methods from the \e Vector class.
 */

#include "vector.h"
#ifdef USE_MPI
#include <mpi.h>
#endif

void Vector::associateChunkData(const int num_elts, int &_halo_start_index,
                                int &_chunk_size, int &_chunk_start_index) {

    _chunk_size = num_elts;
    _chunk_start_index = _halo_start_index;
    _halo_start_index += _chunk_size;
}

void Vector::resize(Dimensions const &in_dims) {

    int imax_loc = in_dims.getNumEltsLoc().i;
    int jmax_loc = in_dims.getNumEltsLoc().j;
    Neighbors ngb_pid = in_dims.getDecomposition().getNgbPid();
    int tmp_halo_start_index;

    dims = in_dims;

    _loc_elts = imax_loc * jmax_loc;
    _halo_elts = countHaloElts(dims);;
    tmp_halo_start_index = _loc_elts;

    rows = _loc_elts + _halo_elts;

    data.resize(rows * cols);

    /*
     * Indentify chunks of halo elements that are stored at the end of the vector
     * and calculate their sizes and starting IndicesBegEnd.
     */
    if (ngb_pid.west != EMPTY)
        associateChunkData(jmax_loc, tmp_halo_start_index,
                           halo_chunk_size.west, halo_chunk_start_index.west);

    if (ngb_pid.south != EMPTY)
        associateChunkData(imax_loc, tmp_halo_start_index,
                           halo_chunk_size.south, halo_chunk_start_index.south);

    if (ngb_pid.north != EMPTY)
        associateChunkData(imax_loc, tmp_halo_start_index,
                           halo_chunk_size.north, halo_chunk_start_index.north);

    if (ngb_pid.east != EMPTY)
        associateChunkData(jmax_loc, tmp_halo_start_index,
                           halo_chunk_size.east, halo_chunk_start_index.east);

    /* Identify which elements should be sent */
    if (ngb_pid.west != EMPTY) {
        on_boarder_ids.west.resize(jmax_loc);
        for(int j = 0; j <= dims.getInternalIndRangeJ().end - dims.getInternalIndRangeJ().beg; ++j) {
            on_boarder_ids.west[j] = j;
        }
    }

    if (ngb_pid.east != EMPTY) {
        on_boarder_ids.east.resize(jmax_loc);
        for(int j = 0; j <= dims.getInternalIndRangeJ().end - dims.getInternalIndRangeJ().beg; ++j) {
            on_boarder_ids.east[j] = j + (dims.getInternalIndRangeJ().end - dims.getInternalIndRangeJ().beg + 1)
                                    * (dims.getInternalIndRangeI().end - dims.getInternalIndRangeI().beg);
        }
    }

    if (ngb_pid.south != EMPTY) {
        on_boarder_ids.south.resize(imax_loc);
        for(int i = 0; i <= dims.getInternalIndRangeI().end - dims.getInternalIndRangeI().beg; ++i) {
            on_boarder_ids.south[i] = (dims.getInternalIndRangeJ().end - dims.getInternalIndRangeJ().beg + 1) * i;
        }
    }

    if (ngb_pid.north != EMPTY) {
        on_boarder_ids.north.resize(imax_loc);
        for(int i = 0; i <= dims.getInternalIndRangeI().end - dims.getInternalIndRangeI().beg; ++i) {
            on_boarder_ids.north[i] = (dims.getInternalIndRangeJ().end - dims.getInternalIndRangeJ().beg)
                                    + (dims.getInternalIndRangeJ().end - dims.getInternalIndRangeJ().beg + 1) * i;
        }
    }
}

void Vector::exchangeRealHalo() {

#ifndef USE_MPI
    // no need to communicate in a non-MPI code
    return;
#else

    int my_rank = getMyRank();
    int num_procs = getNumProcs();
    vector<double> snd_buf_we, rcv_buf_we;                      // Send/receive buffers for the west/east Neighbors
    vector<double> snd_buf_sn, rcv_buf_sn;                      // Send/receive buffers for the south/north Neighbors
    Neighbors ngb_pid = dims.getDecomposition().getNgbPid();
    int imax_loc = dims.getNumEltsLoc().i;
    int jmax_loc = dims.getNumEltsLoc().j;
    MPI_Status status;
    int tag_we = 1;
    int tag_sn = 2;

    // Pre-allocate buffers
    snd_buf_we.resize(jmax_loc);
    rcv_buf_we.resize(jmax_loc);

    snd_buf_sn.resize(imax_loc);
    rcv_buf_sn.resize(imax_loc);

    /* ****************************************************************************************** */
    // Assemble send buffers to west
    if (ngb_pid.west != EMPTY) {
        for(int n = 0; n < halo_chunk_size.west; ++n) {
            snd_buf_we[n] = data[on_boarder_ids.west[n]];
        }
        MPI_Send(snd_buf_we.data(), snd_buf_we.size(), MPI_DOUBLE, ngb_pid.west, tag_we, MPI_COMM_WORLD);
    }

    // Assemble send buffers to south
    if (ngb_pid.south != EMPTY) {
        for(int n = 0; n < halo_chunk_size.south; ++n) {
            snd_buf_sn[n] = data[on_boarder_ids.south[n]];
        }
        MPI_Send(snd_buf_sn.data(), snd_buf_sn.size(), MPI_DOUBLE, ngb_pid.south, tag_sn, MPI_COMM_WORLD);
    }

    // Receive from east
    if (ngb_pid.east != EMPTY) {
        MPI_Recv(rcv_buf_we.data(), rcv_buf_we.size(), MPI_DOUBLE, ngb_pid.east, tag_we, MPI_COMM_WORLD, &status);
        int id = halo_chunk_start_index.east;
        for(int n = 0; n < halo_chunk_size.east; ++n) {
            data[id + n] = rcv_buf_we[n];
        }
    }

    // Receive from north
    if (ngb_pid.north != EMPTY) {
        MPI_Recv(rcv_buf_sn.data(), rcv_buf_sn.size(), MPI_DOUBLE, ngb_pid.north, tag_sn, MPI_COMM_WORLD, &status);
        int id = halo_chunk_start_index.north;
        for(int n = 0; n < halo_chunk_size.north; ++n) {
            data[id + n] = rcv_buf_sn[n];
        }
    }
    /* ****************************************************************************************** */
    /* ****************************************************************************************** */
    // Assemble send buffers to east
    if (ngb_pid.east != EMPTY) {
        for(int n = 0; n < halo_chunk_size.east; ++n) {
            snd_buf_we[n] = data[on_boarder_ids.east[n]];
        }
        MPI_Send(snd_buf_we.data(), snd_buf_we.size(), MPI_DOUBLE, ngb_pid.east, tag_we, MPI_COMM_WORLD);
    }

    // Assemble send buffers to north
    if (ngb_pid.north != EMPTY) {
        for(int n = 0; n < halo_chunk_size.north; ++n) {
            snd_buf_sn[n] = data[on_boarder_ids.north[n]];
        }
        MPI_Send(snd_buf_sn.data(), snd_buf_sn.size(), MPI_DOUBLE, ngb_pid.north, tag_sn, MPI_COMM_WORLD);
    }

    // Receive from west
    if (ngb_pid.west != EMPTY) {
        MPI_Recv(rcv_buf_we.data(), rcv_buf_we.size(), MPI_DOUBLE, ngb_pid.west, tag_we, MPI_COMM_WORLD, &status);
        int id = halo_chunk_start_index.west;
        for(int n = 0; n < halo_chunk_size.west; ++n) {
            data[id + n] = rcv_buf_we[n];
        }
    }

    // Receive from south
    if (ngb_pid.south != EMPTY) {
        MPI_Recv(rcv_buf_sn.data(), rcv_buf_sn.size(), MPI_DOUBLE, ngb_pid.south, tag_sn, MPI_COMM_WORLD, &status);
        int id = halo_chunk_start_index.south;
        for(int n = 0; n < halo_chunk_size.south; ++n) {
            data[id + n] = rcv_buf_sn[n];
        }
    }
    /* ****************************************************************************************** */

    MPI_Barrier(MPI_COMM_WORLD);
#endif
}
