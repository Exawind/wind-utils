//  Copyright 2016 National Renewable Energy Laboratory
//
//  Licensed under the Apache License, Version 2.0 (the "License");
//  you may not use this file except in compliance with the License.
//  You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
//

#ifndef PARALLELINFO_H
#define PARALLELINFO_H

#include <string>
#include <iostream>

#include "mpi.h"

namespace sierra {
namespace nalu {

namespace wind_utils {

class NullBuffer : public ::std::streambuf
{
public:
    int overflow(int c) { return c; }
};

/** Utility class to track MPI information
 */
class ParallelInfo
{
public:
    ParallelInfo()
        : comm_(MPI_COMM_WORLD),
          null_stream_(&null_buffer_)
    {
        MPI_Comm_size(comm_, &size_);
        MPI_Comm_rank(comm_, &rank_);
    }

    ~ParallelInfo() = default;

    //! Communicator instance
    MPI_Comm comm() const { return comm_; }

    //! Rank of the running process
    int rank() const { return rank_ ; }

    //! Total size of the MPI communicator
    int size() const { return size_ ; }

    //! Flag indicating whether this is rank 0
    bool master() const { return (rank_ == 0); }

    //! Output stream that only prints on the master rank
    std::ostream& info() const
    { return (rank_ == 0) ? std::cout : null_stream_; }

    /** Output stream that prints on every MPI rank
     *
     *  If prefix is true, then it will print the MPI rank so that the messages
     *  are labeled by the processor rank.
     */
    std::ostream& pInfo(bool prefix = true) const
    {
        if (prefix)
            std::cout << std::endl << "[" << rank_ << "] ";
        return std::cout;
    }

private:
    //! Parallel communicator instance
    MPI_Comm comm_;

    //! Current processor rank
    int rank_{0};

    //! Total number of MPI ranks
    int size_{1};

    mutable NullBuffer null_buffer_;

    mutable std::ostream null_stream_;
};

}  // wind_utils

/** Get the global MPI info object
 */
inline const wind_utils::ParallelInfo& get_mpi()
{
    static wind_utils::ParallelInfo pinfo;
    return pinfo;
}

}  // nalu
}  // sierra


#endif /* PARALLELINFO_H */
