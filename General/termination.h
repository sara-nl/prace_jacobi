//
// Created by Maxim Masterov on 09/06/2021.
//

#ifndef TERMINATION_H
#define TERMINATION_H

#ifdef USE_MPI
#include <mpi.h>
#endif

/*!
 * @brief Terminate the program.
 */
inline void terminateExecution() {
#ifdef USE_MPI
    MPI_Abort(MPI_COMM_WORLD, 1);
#else
    exit(1);
#endif
}

#endif //TERMINATION_H
