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
 * @file solver.cpp
 * @brief Contains definitions of methods from the \e Solver class.
 */

#include "solver.h"

void Solver::copyVector(Vector &vec_in, Vector &vec_out) {

    /*
     * for every n-th elements in `vec_in`
     *     assign element of vec_in(n) to vec_out(n)
     */
    // NOT_IMPLEMENTED

    for(int i = 0; i < vec_in.numRows(); ++i) {
        vec_out(i) = vec_in(i);
    }    
}

void Solver::calculateResidual(Matrix &A, Vector &x, Vector &b, Vector &res) {

    /*
     * assign `b` to `res`
     *
     * for vector `x` and matrix `A`, the residual `res`
     * is calculated as:
     *     res(i) = b(i) - sum( A(i, j) * x(j) )
     */
    // NOT_IMPLEMENTED

    copyVector(b, res);

    for(int i = 0; i < A.numRows(); ++i) {
        double sum = 0.0;
        for(int j = 0; j < A.numCols(); ++j) {
            sum += A(i, j) * x(j);
        }
        res(i) -= sum;
    }
}

double Solver::calculateNorm(Vector &vec) {

    /*
     * for vector `vec` with n elements
     *   L2-norm = sqrt( sum( vec(n) * vec(n) ) )
     *             |     \_____________________/|
     *             | local for every MPI process|
     *              \__________________________/
     *            one value for all MPI processes
     */
    // NOT_IMPLEMENTED
    
    double sum = 0.0;

    for(int i = 0; i < vec.getLocElts(); ++i) {
        sum += vec(i) * vec(i);
    }

    findGlobalSum(sum);

    return sqrt(sum);
}

void Solver::solveJacobi(Matrix &A, Vector &x, Vector &b) {

    int iter = 0;                   // Iteration counter
    int max_iter = 10000;           // Maximum number of iterations
    double tolerance = 1e-6;        // Stopping criteria
    double omega = 2./3.;            // Under-relaxation factor
    double residual_norm = 0.0;     // Normalized residual
    Vector x_old;                   // Old solution
    Vector res;                     // Residual vector
    int my_rank = 0;                // Process rank (0 in non-MPI case)

    my_rank = getMyRank();

    x_old.resize(x.getDimensions());
    res.resize(x.getDimensions());

    residual_norm = 10. * tolerance;

    copyVector(x, x_old);

    /* Start the main loop */
#ifdef USE_MPI
    /*
     * Note, that the data should be exchanged between the real and halo cells.
     * Place the call for `x.exchangeRealHalo();` at the right place in the
     * iterative loop.
     */
    // NOT_IMPLEMENTED
#endif
    while ( (iter < max_iter) && (residual_norm > tolerance) ) {
        
        for(int i = A.numRows() - 1; i >= 0; i--) {
            double diag = 1.;          // Diagonal element
            double sigma = 0.0;        // Just a temporary value

            x(i) = b(i);

            for(int j = 0; j < A.numCols(); ++j) {
                if (j != i)
                    sigma = sigma + A(i, j) * x_old(j);
                else
                    diag = A(i, j);
            }
            x(i) = (x(i) - sigma) * omega / diag;
        }
        
        x.exchangeRealHalo();

        for(int i = 0; i < x.numRows(); ++i) {
            x(i) += (1 - omega) * x_old(i);
            x_old(i) = x(i);
        }

        calculateResidual(A, x, b, res);
        residual_norm = calculateNorm(res) / calculateNorm(b);

        if (my_rank == 0)
            cout << iter << '\t' << residual_norm << endl;

        ++iter;
    }
}
