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

#include "PerfUtils.h"
#include <sstream>

namespace sierra {
namespace nalu {

std::string human_bytes(const size_t bytes)
{
    std::ostringstream hbytes;
    const double kb = 1024.0;
    const double mb = kb * kb;
    const double gb = mb * kb;

    const double tmp = static_cast<double>(bytes);

    if (tmp >= gb)
        hbytes << std::setw(8) << std::fixed << std::setprecision(3)
               << tmp / gb << " GB";
    else if (tmp >= mb)
        hbytes << std::setw(8) << std::fixed << std::setprecision(3)
               << tmp / mb << " MB";
    else if (tmp >= kb)
        hbytes << std::setw(8) << std::fixed << std::setprecision(3)
               << tmp / kb << " KB";
    else
        hbytes << std::setw(8) << std::fixed << std::setprecision(3)
               << tmp << "  B";

    return hbytes.str();
}

void summarize_memory_usage(MPI_Comm comm, std::ostream& out)
{
    int rank = stk::parallel_machine_rank(comm);
    size_t hwm_max =0, hwm_min = 0, hwm_avg = 0;

    stk::get_memory_high_water_mark_across_processors(
        comm, hwm_max, hwm_min, hwm_avg);

    if (rank == 0)
        out << "\nMemory usage: "
            << "Avg: " << human_bytes(hwm_avg)
            << "; Min: " << human_bytes(hwm_min)
            << "; Max: " << human_bytes(hwm_max)
            << std::endl << std::endl;
}

}  // nalu
}  // sierra
