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
 * @file field.cpp
 * @brief Contains definition of methods from the \e Field class.
 */

#include "field.h"

void Field::resize(Dimensions &in_dims) {

    int imax_loc = in_dims.getNumEltsLoc().i;
    int jmax_loc = in_dims.getNumEltsLoc().j;

    dims = in_dims;

    _loc_elts = imax_loc * jmax_loc;
    _halo_elts = countHaloElts(in_dims);

    rows = imax_loc;
    cols = jmax_loc;

    data.resize(rows * cols);

    enumerateIDs();
}

void Field::enumerateIDs() {

    int counter = 0;
    int total_num_elts = 0;

    total_num_elts = rows * cols + _halo_elts;

    /*
     * Check corners and add them to the element counter. In general 2d case the
     * corners may appear at the following four sites:
     *  - west-south
     *  - west-north
     *  - east-south
     *  - east-north
     */
    if (dims.getDecomposition().getNgbPid().west != EMPTY &&
            dims.getDecomposition().getNgbPid().south != EMPTY) {
        ++total_num_elts;
    }
    if (dims.getDecomposition().getNgbPid().west != EMPTY &&
            dims.getDecomposition().getNgbPid().north != EMPTY) {
        ++total_num_elts;
    }
    if (dims.getDecomposition().getNgbPid().east != EMPTY &&
            dims.getDecomposition().getNgbPid().south != EMPTY) {
        ++total_num_elts;
    }
    if (dims.getDecomposition().getNgbPid().east != EMPTY &&
            dims.getDecomposition().getNgbPid().north != EMPTY) {
        ++total_num_elts;
    }

    ids.resize(total_num_elts);

    /*
     * Initialize vactor of ids with "EMPTY". This will help to indicate corner
     * halo elements
     */
    for(int i = 0; i < total_num_elts; ++i) {
        ids[i] = EMPTY;
    }

    /* Enumerate internal elements */
    for(int i = dims.getInternalIndRangeI().beg; i <= dims.getInternalIndRangeI().end; ++i) {
        for(int j = dims.getInternalIndRangeJ().beg; j <= dims.getInternalIndRangeJ().end; ++j) {
            ids[j + i * dims.getNumElts().j] = counter;
            ++counter;
        }
    }

    /* Enumerate halo elements in the west */
    if (dims.getDecomposition().getNgbPid().west != EMPTY) {
        for(int j = dims.getInternalIndRangeJ().beg; j <= dims.getInternalIndRangeJ().end; ++j) {
            ids[j] = counter;
            ++counter;
        }
    }

    /* Enumerate halo elements in the south */
    if (dims.getDecomposition().getNgbPid().south != EMPTY) {
        for(int i = dims.getInternalIndRangeI().beg; i <= dims.getInternalIndRangeI().end; ++i) {
            ids[i * dims.getNumElts().j] = counter;
            ++counter;
        }
    }

    /* Enumerate halo elements in the north */
    if (dims.getDecomposition().getNgbPid().north != EMPTY) {
        for(int i = dims.getInternalIndRangeI().beg; i <= dims.getInternalIndRangeI().end; ++i) {
            ids[dims.getNumElts().j - 1 + i * dims.getNumElts().j] = counter;
            ++counter;
        }
    }

    /* Enumerate halo elements in the east */
    if (dims.getDecomposition().getNgbPid().east != EMPTY) {
        for(int j = dims.getInternalIndRangeJ().beg; j <= dims.getInternalIndRangeJ().end; ++j) {
            ids[j + (dims.getNumElts().i - 1) * dims.getNumElts().j] = counter;
            ++counter;
        }
    }
}
