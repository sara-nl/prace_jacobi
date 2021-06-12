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

#ifdef USE_GPU
#pragma omp declare target
void Solver::copyVector(double *vec_in, double *vec_out, int num_rows) {

    /*
     * for every n-th elements in `vec_in`
     *     assign element of vec_in(n) to vec_out(n)
     */
    // NOT_IMPLEMENTED
#pragma omp distribute parallel for simd
    for(int i = 0; i < num_rows; ++i) {
        vec_out[i] = vec_in[i];
    }    
}
#pragma omp end declare target

#pragma omp declare target
void Solver::calculateResidual(double *A, double *x, double *b, double *res,
                               int num_rows, int num_cols) {

    /*
     * assign `b` to `res`
     *
     * for vector `x` and matrix `A`, the residual `res`
     * is calculated as:
     *     res(i) = b(i) - sum( A(i, j) * x(j) )
     */
    // NOT_IMPLEMENTED

    copyVector(b, res, num_rows);

#pragma omp distribute parallel for
    for(int i = 0; i < num_rows; ++i) {
        double sum = 0.0;
#pragma omp simd reduction(+:sum)
        for(int j = 0; j < num_cols; ++j) {
            sum += A[j + i * num_cols] * x[j];
        }
        res[i] -= sum;
    }
}
#pragma omp end declare target

#pragma omp declare target
double Solver::calculateNorm(double *vec, int num_loc_elts) {

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

#pragma omp parallel for simd reduction(+ : sum)
    for(int i = 0; i < num_loc_elts; ++i) {
        sum += vec[i] * vec[i];
    }

    // findGlobalSum(sum);

    return sqrt(sum);
}
#pragma omp end declare target

void Solver::solveJacobiGPU(Matrix &A_obj, Vector &x_obj, Vector &b_obj) {

    int iter = 0;                   // Iteration counter
    int max_iter = 10000;           // Maximum number of iterations
    double tolerance = 1e-6;        // Stopping criteria
    double omega = 2./3.;            // Under-relaxation factor
    double residual_norm = 0.0;     // Normalized residual
    Vector x_old_obj;                   // Old solution
    Vector res_obj;                     // Residual vector
    int my_rank = 0;                // Process rank (0 in non-MPI case)
    double norm_b;
    
    my_rank = getMyRank();

    x_old_obj.resize(x_obj.getDimensions());
    res_obj.resize(x_obj.getDimensions());

    residual_norm = 10. * tolerance;

    // Define raw pointer that will be used in the target region
    double *x = x_obj.getData();
    double *x_old = x_old_obj.getData();
    double *b = b_obj.getData();
    double *res = res_obj.getData();
    double *A = A_obj.getData();
    int num_rows = A_obj.numRows();
    int num_cols = A_obj.numCols();
    int num_elts_A = num_rows * num_cols;
    int num_loc_elts_x = x_obj.numRows();

    copyVector(x, x_old, num_rows);
    
    /* Start the main loop */
#ifdef USE_MPI
    /*
     * Note, that the data should be exchanged between the real and halo cells.
     * Place the call for `x.exchangeRealHalo();` at the right place in the
     * iterative loop.
     */
    // NOT_IMPLEMENTED
#endif

#pragma omp target data map(to: A[0:num_elts_A], b[0:num_loc_elts_x], res[0:num_loc_elts_x]) \
                        map(tofrom: x[0:num_loc_elts_x], x_old[0:num_loc_elts_x])
{
#pragma omp target map(from: norm_b)
    norm_b = calculateNorm(b, num_loc_elts_x);

    while ( (iter < max_iter) && (residual_norm > tolerance) ) {

#pragma omp target teams
#pragma omp distribute parallel for
        for(int i = num_rows - 1; i >= 0; i--) {
            double diag = 1.;          // Diagonal element
            double sigma = 0.0;        // Just a temporary value

            x[i] = b[i];
          
#pragma omp simd reduction(+ : sigma)
            for(int j = 0; j < num_cols; ++j) {
                sigma += A[j + i * num_cols] * x_old[j];
            }
            diag = A[i + i * num_cols];
            sigma -= diag * x_old[i];
            x[i] = (x[i] - sigma) * omega / diag;
        }

        // x.exchangeRealHalo();

#pragma omp target teams
#pragma omp distribute parallel for simd
        for(int i = 0; i < num_loc_elts_x; ++i) {
            x[i] += (1 - omega) * x_old[i];
            x_old[i] = x[i];
        }

#pragma omp target teams
        calculateResidual(A, x, b, res, num_cols, num_rows);

#pragma omp target map(tofrom:residual_norm)
        residual_norm = calculateNorm(res, num_loc_elts_x);

        residual_norm = residual_norm / norm_b;

        if (my_rank == 0)
            cout << iter << '\t' << residual_norm << endl;

        ++iter;
    }
}

// Download the solution back from the device
#pragma omp target update from(x[0:num_loc_elts_x])

}
#endif