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

#ifdef USE_MPI
void Vector::sendAndReceive(vector<double> &snd_buf, vector<double> &rcv_buf, 
                    int ngb_pid, int tag, MPI_Request *request) {

    MPI_Isend(snd_buf.data(), snd_buf.size(), MPI_DOUBLE, ngb_pid, tag, MPI_COMM_WORLD, &request[0]);
    MPI_Irecv(rcv_buf.data(), rcv_buf.size(), MPI_DOUBLE, ngb_pid, tag, MPI_COMM_WORLD, &request[1]);
}

void Vector::packSndBuffer(vector<double> &snd_buf, std::vector<int> &on_boarder_ids, int halo_chunk_size) {
    
    for(int n = 0; n < halo_chunk_size; ++n) {
        snd_buf[n] = data[on_boarder_ids[n]];
    }
}

void Vector::unpackRcvBuffer(vector<double> &rcv_buf, int halo_chunk_start_index, 
                             int halo_chunk_size, int ngb_pid) {
    
    if (ngb_pid != EMPTY) {
        for(int n = 0; n < halo_chunk_size; ++n) {
            data[halo_chunk_start_index + n] = rcv_buf[n];
        }
    }
}

void Vector::communicateToNgb(vector<double> &snd_buf, vector<double> &rcv_buf,
                              vector<int> &on_boarder_ids, int halo_chunk_size,
                              int ngb_pid, int tag, MPI_Request *request) {
    if (ngb_pid != EMPTY) {
        packSndBuffer(snd_buf, on_boarder_ids, halo_chunk_size);
        sendAndReceive(snd_buf, rcv_buf, ngb_pid, tag, request);
    }
}

void Vector::waitForAll(Neighbors &ngb_pid, 
                MPI_Request *request_w, MPI_Request *request_e,
                MPI_Request *request_s, MPI_Request *request_n) {
    
    MPI_Status snd_status[2];
    if (ngb_pid.west != EMPTY) {
        MPI_Waitall(2, request_w, snd_status);
    }
    if (ngb_pid.east != EMPTY) {
        MPI_Waitall(2, request_e, snd_status);
    }
    if (ngb_pid.south != EMPTY) {
        MPI_Waitall(2, request_s, snd_status);
    }
    if (ngb_pid.north != EMPTY) {
        MPI_Waitall(2, request_n, snd_status);
    }
}
#endif

//
// Using MPI_Isend() and MPI_Irecv()
//
void Vector::exchangeRealHalo() {

#ifndef USE_MPI
    // no need to communicate in a non-MPI code
    return;
#else
    vector<double> snd_buf_w, rcv_buf_w;                      // Send/receive buffers for the west neighbors
    vector<double> snd_buf_e, rcv_buf_e;                      // Send/receive buffers for the east neighbors
    vector<double> snd_buf_s, rcv_buf_s;                      // Send/receive buffers for the south neighbors
    vector<double> snd_buf_n, rcv_buf_n;                      // Send/receive buffers for the north neighbors

    Neighbors ngb_pid = dims.getDecomposition().getNgbPid();
    int imax_loc = dims.getNumEltsLoc().i;
    int jmax_loc = dims.getNumEltsLoc().j;
    MPI_Status status;
    int tag_we = 1;
    int tag_sn = 2;
    int num_ngb = 4;
    MPI_Request request_w[2], request_e[2], request_s[2], request_n[2];

    // Pre-allocate buffers
    snd_buf_w.resize(jmax_loc);
    snd_buf_e.resize(jmax_loc);
    rcv_buf_w.resize(jmax_loc);
    rcv_buf_e.resize(jmax_loc);

    snd_buf_s.resize(imax_loc);
    snd_buf_n.resize(imax_loc);
    rcv_buf_s.resize(imax_loc);
    rcv_buf_n.resize(imax_loc);

    // Communicate with the neighbors on all sides
    communicateToNgb(snd_buf_w, rcv_buf_w, on_boarder_ids.west, halo_chunk_size.west, 
                     ngb_pid.west, tag_we, request_w);
    communicateToNgb(snd_buf_e, rcv_buf_e, on_boarder_ids.east, halo_chunk_size.east,
                     ngb_pid.east, tag_we, request_e);
    communicateToNgb(snd_buf_s, rcv_buf_s, on_boarder_ids.south, halo_chunk_size.south,
                     ngb_pid.south, tag_sn, request_s);
    communicateToNgb(snd_buf_n, rcv_buf_n, on_boarder_ids.north, halo_chunk_size.north,
                     ngb_pid.north, tag_sn, request_n);

    // Wait for all communications to finish
    waitForAll(ngb_pid, request_w, request_e, request_s, request_n);
    
    // Unpack received buffers
    unpackRcvBuffer(rcv_buf_w, halo_chunk_start_index.west, halo_chunk_size.west, ngb_pid.west);
    unpackRcvBuffer(rcv_buf_e, halo_chunk_start_index.east, halo_chunk_size.east, ngb_pid.east);
    unpackRcvBuffer(rcv_buf_s, halo_chunk_start_index.south, halo_chunk_size.south, ngb_pid.south);
    unpackRcvBuffer(rcv_buf_n, halo_chunk_start_index.north, halo_chunk_size.north, ngb_pid.north);

#endif
}
