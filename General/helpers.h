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
 * @file helpers.h
 * @brief Contains declarations of miscellaneous functions.
 * These functions are not essential and are used to simplify the main code.
 */

#ifndef HELPERS_H_
#define HELPERS_H_

#include "dimensions.h"

class Helpers {
public:
    /*!
     * @brief Evaluate parameters passed through the command line.
     *
     * Evaluates parameters and assignes them to the number of cells in the domain.
     *
     * Note, if one parameter is passed it will be used for both directions.
     *
     * @param argc [in] Number of command line arguments
     * @param argv [in] Vector of command line arguments
     * @param dims [out] Structure with Dimensions of the domain
     */
    void setDimensionsAndDecompose(int argc, char** argv, Dimensions &dims);

    /*!
     * @brief Start the timer and return the current time (in seconds) starting
     *        from "some event in the past".
     */
    double tic();

    /*!
     * @brief Stop the timer and return the current time (in seconds) starting
     *        from "some event in the past".
     */
    double toc();

    /*!
     * @brief Parse the input parameters from CL.
     * @param argc Number of CL parameters.
     * @param argv CL parameters.
     * @param elts_glob Number of global elements in each direction.
     * @param num_procs Number of local elements in each direction.
     */
    void parseInput(int argc, char** argv, IndicesIJ &elts_glob, IndicesIJ &num_procs);

private:
    /*!
     * @brief Terminate execution due to error in the input parameters.
     */
    void terminateDueToParserFailure();
};

#endif /* HELPERS_H_ */
