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

#ifndef STRUCTS_H
#define STRUCTS_H

/*!
 * @brief Structure of begin/end indices.
 */
struct IndicesBegEnd {
    int beg = 0;
    int end = 0;

    IndicesBegEnd() { }
    IndicesBegEnd(int _beg, int _end) : beg(_beg), end(_end) { }
};

/*!
 * @brief Structure of {i,j} indices.
 */
struct IndicesIJ {
    int i = 0;
    int j = 0;

    IndicesIJ() { }
    IndicesIJ(int _i, int _j) : i(_i), j(_j) { }
};

/*!
 * @brief Structure of face values.
 * Used to set boundary conditions or interpolation weight.
 */
struct Faces {
    double east = 0.0;
    double west = 0.0;
    double south = 0.0;
    double north = 0.0;
    double central = 0.0;
};

/*!
 * @brief Structure of the neighboring processes ids.
 * @note \e EMPTY means there is no neighbor in that direction.
 */
struct Neighbors {
    int east = EMPTY;
    int west = EMPTY;
    int south = EMPTY;
    int north = EMPTY;
    int central = EMPTY;
};
#endif
