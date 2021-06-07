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
 * @file common.h
 * @brief Contains C-based declarations of the main MPI functions
 */

#ifndef COMMON_H
#define COMMON_H

#include <string>

#ifdef USE_MPI
#include <mpi.h>
#endif

/*!
 * @brief Find minimum value among all processes.
 * @note This function replaces input with the output.
 * @param value [in/out] The value to find.
 */
void findGlobalMin(double &value);

/*!
 * @brief Find maximum value among all processes.
 * @note This function replaces input with the output.
 * @param value [in/out] The value to find.
 */
void findGlobalMax(double &value);

/*!
 * @brief Perform global summation of all elements and returns a single (global)
 *        result.
 * @note This function replaces input with the output.
 * @param value [in/out] The value to sum.
 */
void findGlobalSum(double &value);
void findGlobalSum(int &value);

/*!
 * @brief Get the rank of the local process.
 */
int getMyRank();

/*!
 * @brief Get the number of processes.
 */
int getNumProcs();

/*!
 * @brief Initialize MPI calls.
 * @param argc [in] Number of command line arguments.
 * @param argv [in] List of command line arguments.
 */
void initialize(int argc, char** argv);

/*!
 * @brief Finalize MPI calls.
 */
void finalize();

/*!
 * @brief Terminate the program.
 */
void terminateExecution();

/*!
 * @brief Print the message to the terminal.
 * @param str [in] The message.
 */
void printByRoot(const std::string& str);
#endif
