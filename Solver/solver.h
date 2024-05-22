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
 * @file solver.h
 * @brief Contains declaration of the \e Solver class.
 */

#ifndef SOLVER_H_
#define SOLVER_H_

#ifdef USE_MPI
#include <mpi.h>
#endif
#include <vector>
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>

#include "../DataTypes/matrix.h"
#include "../DataTypes/vector.h"
#include "../General/structs.h"

using namespace std;

/*!
 * @class Solver
 * @brief Responsible for all math operations.
 */
class Solver {
public:
    /*!
     * @brief Calculate the residual \f[ r = b - Ax \f].
     * @param A [in] Matrix
     * @param x [in] Vector of unknowns
     * @param b [in] Vector of right hand side
     * @param res [out] Vector of residual
     */
    void calculateResidual(Matrix &A, Vector &x, Vector &b, Vector &res);

    /*!
     * @brief Calculate the L2-norm.
     * @param vec [in] Vector
     * @return Value of L2-norm
     */
    double calculateNorm(Vector &vec);

    /*!
     * @brief Copy elements of one vector to another vector.
     * @param [in] vec_in Vector to be copied from
     * @param [out] vec_out Vector to be copied to
     */
    void copyVector(Vector &vec_in, Vector &vec_out);

    /*!
     * @brief Solve the provided linear system \f[ A x = b \f] using Jacobi solver.
     * @note Memory for the vectors and matrix should be pre-allocated.
     * @param A [in] Matrix
     * @param x [out] Vector of unknowns
     * @param b [in] Vector of right hand side
     */
    void solveJacobi(Matrix &A, Vector &x, Vector &b);
};



#endif /* SOLVER_H_ */
