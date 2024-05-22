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
 * @file system.h
 * @brief Contains declaration of the System class
 */

#ifndef SYSTEM_H_
#define SYSTEM_H_

#include "../DataTypes/matrix.h"
#include "../DataTypes/vector.h"
#include "../DataTypes/field.h"
#include "../General/structs.h"

/*!
 * @class System
 * @brief Responsible for setup of the numerical problem
 */
class System {

public:
    /*!
     * @brief Allocate memory for the linear system and the fied.
     * @param dims [in] Structure with Dimensions of the domain
     * @param T [out] Field of temperature
     * @param A [out] Matrix
     * @param x [out] Vector of unknowns
     * @param b [out] Vector of right hand side
     */
    void allocateMemory(Dimensions &dims, Field &T,
                        Matrix &A, Vector &x, Vector &b);

    /*!
     * @brief Assemble the linear system of a form \f[ A x = b \f].
     *
     * @note The system will be assembled with Dirichlet boundary conditions at all walls.
     *
     * @param bondary_values [in] Structure with boundary values
     * @param T [in] Field
     * @param A [out] Matrix
     * @param x [out] Vector of unknowns
     * @param b [out] Vector of right hand side
     */
    void assembleSystem(Faces &bondary_values, Field &T,
                        Matrix &A, Vector &x, Vector &b);

    /*!
     * @brief Copy the solution of the linear system back to the field.
     * @param x [in] Vector of unknowns
     * @param T [out] Field of temperature
     */
    void copySolution(Vector &x, Field &T);
};

#endif /* SYSTEM_H_ */
