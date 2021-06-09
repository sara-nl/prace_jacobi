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
 * @file matrix.h
 * @brief Contains declaration of the \e Matrix class.
 */

#ifndef MATRIX_H
#define MATRIX_H

#ifdef USE_MPI
#include <mpi.h>
#endif
#include <iostream>
#include <vector>
#include "../General/dimensions.h"

using namespace std;

/*!
 * @class Matrix
 * @brief Represents dense matrix.
 * The matrix can be local or distributed.
 */
class Matrix {
protected:
    vector<double> data;    // Vector of elements
    int rows;               // Number of rows
    int cols;               // Number of columns

    int _loc_elts;          // Number of local elements
    int _halo_elts;         // Number of halo elements

    Dimensions dims;        // Dimensions of the numerical domain

public:
    /*!
     * @brief Default constructor.
     */
    Matrix() : rows(0), cols(0), _loc_elts(0), _halo_elts(0) {}

    /*!
     * @brief Default destructor.
     */
    virtual ~Matrix() { }

    /*!
     * @brief Allocate memory for the matrix.
     * @note The memory is allocated but not initialized.
     * @param in_dims [in] Dimensions of the numerical problem.
     */
    virtual void resize(Dimensions const &in_dims);

    /*!
     * @brief Print the matrix.
     * If matrix is distributer, the function prints it sequentially process
     * by process.
     */
    void print();

    /*!
     * @brief Return a reference to the element at the specified IndicesBegEnd: row
     *        and column.
     * @param row [in] Row.
     * @param col [in] Column.
     */
    inline double &operator()(int row, int col) {
        return data[col + row * cols];
    }

    /*!
     * @brief Return raw data.
     */
    inline double *getData() {
        return data.data();
    }

    /*!
     * @brief Return number of local rows including rows which represent halo 
     *        elements, if there are any.
     * @note If this method is called for the object of inherited \e Field class,
     * the method will return the number of elements in i-th direction.
     */
    inline int numRows() {
        return rows;
    }

    /*!
     * @brief Return the total number of local columns including columns which
     *        represent halo elements.
     * @note If this method is called for the object of inherited \e Field class,
     * the method will return the number of elements in j-th direction.
     */
    inline int numCols() {
        return cols;
    }

    /*!
     * @brief Return the total number of local elements in the matrix including
     *        halo elements.
     */
    inline int size() {
        return rows * cols;
    }

    /*!
     * @brief Return a copy of a structure of Dimensions.
     */
    inline const Dimensions &getDimensions() const {
        return dims;
    }

    /**********************************/
    /********** Just helpers **********/
    /**********************************/
    /*!
     * @brief Return the number of local elements excluding halo elements.
     */
    inline int getLocElts() {
        return _loc_elts;
    }

    /*!
     * @brief Return the number of local halo elements.
     */
    inline int getHaloElts() {
        return _halo_elts;
    }

protected:
    /*!
     * @brief Check for existence of the neighboring processes and counts
     *        the number of halo cells.
     * @param in_dims [in] Structure of Dimensions of the numerical domain.
     * @return The total number of halo elements.
     */
    int countHaloElts(Dimensions const &in_dims);
};

#endif
