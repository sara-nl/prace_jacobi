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

/*!
 * @file field.h
 * @brief Contains declaration of the \e Field class.
 */

#ifndef FIELD_H
#define FIELD_H

#include "matrix.h"
#include "../General/structs.h"

using namespace std;

/*!
 * @class Field
 * @brief Represents generic 2d field.
 * The field can be local or distributed.
 */
class Field : public Matrix {

    vector<int> ids;                // Vector of IDs: first enumerates internal
                                    // cells, than halo cells, than corner cells

public:
    /*!
     * @brief Default constructor
     */
    Field() { }

    /*!
     * @brief Allocate memory.
     * @note The allocated memory is not initialized.
     * @param in_dims [in] Dimensions of the numerical problem.
     */
    void resize(Dimensions &in_dims);

    /*!
     * @brief Return a vector of IDs.
     */
    inline const vector<int> &getIDs() const {
        return ids;
    }

    /*!
     * @brief Return ID of element from its IndicesBegEnd.
     * @param i Index in i-th direction.
     * @param j Index in j-th direction.
     */
    inline int &getID(const int i, const int j) {
        return ids[j + i * dims.getNumElts().j];
    }

private:
    /*!
     * @brief Fill the vector of IDs.
     */
    void enumerateIDs();
};

#endif
