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
 * @file io.h
 * @brief Contains declaration of the \e IO class.
 */

#ifndef IO_IO_H_
#define IO_IO_H_

#ifdef USE_MPI
#include <mpi.h>
#endif
#include <vector>
#include <string>

#include "../DataTypes/field.h"
#include "../General/dimensions.h"

/*!
 * @class IO
 * @brief Responsible for IO operations with the file system.
 */
class IO {

public:
    /*!
     * @brief Write data into the file.
     * @param dims [in] Structure with Dimensions of the domain
     * @param T [in] Field of temperature
     */
    void writeFile(std::string file_name, Dimensions &dims, Field &T);

private:
#ifdef USE_MPI
    /*!
     * @brief Write data into the file by a single process.
     * @param mpi_file [in] MPI file handler
     * @param dims [in] Structure with Dimensions of the domain
     * @param T [in] Field of temperature
     */
    void writeByRoot(MPI_File &mpi_file, Dimensions &dims, Field &T);

    /*!
     * @brief Write data into the file by all process in parallel.
     * @param mpi_file [in] MPI file handler
     * @param dims [in] Structure with Dimensions of the domain
     * @param T [in] Field of temperature
     */
    void writeByAll(MPI_File &mpi_file, Dimensions &dims, Field &T);

    /*!
     * @brief Generate a numerical grid and store result in a form of 1D array.
     * @param dims [in] Structure with Dimensions of the domain
     * @param T [in] Field of temperature
     * @param grid_1D [out] The generated grid
     */
    void generateGrid(Dimensions &dims, Field &T, vector<double> &grid_1D);

    /*!
     * @brief Convert a 2D field to a 1D field.
     * @param field_2D [in] 2D field
     * @param field_1D [in] 2D field
     */
    void convertTo1D(Field &field_2D, vector<double> &field_1D);

    /*!
     * @brief Garther a field by the root process.
     * @param root_pid [in] PID of the root process
     * @param field [in] Field to be gathered
     * @param rcv_buffer [out] Received buffer that contains the gathered field
     */
    void gatherField(int root_pid, vector<double> &field, vector<double> &rcv_buffer);
#endif
};



#endif /* IO_IO_H_ */
