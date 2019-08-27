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

#include <iostream>

#include "struct_grid/StructGrid.h"
#include "struct_grid/StructGridIx.h"
#include "core/ParallelInfo.h"
#include "stk_util/parallel/Parallel.hpp"

int main(int argc, char** argv)
{
    namespace sgix = sierra::nalu::sgix;
    const auto comm = stk::parallel_machine_init(&argc, &argv);

    sierra::nalu::StructGrid grid;

    grid.set_num_ghost(1);
    grid.set_global_grid(64, 32, 32);
    grid.set_partitions(sierra::nalu::get_mpi().size(), 1, 1);

    // const auto& global = grid.global();

    // const sgix::Indexer<sgix::LeftLayout> ixl(global);
    // const sgix::Indexer<sgix::RightLayout> ixr(global);

    // std::cout << ixl(1, 0, 1) << " " << ixr(0, 1, 1) << std::endl;

    // const sgix::PeriodicIndexer<sgix::LeftLayout> ixl(global);

    // std::cout << ixl(0, 0, -1) << " " << ixl(64, 1, 0) << std::endl;

    const auto& local = grid.local();
    std::cout << local.start[0] << " " << local.start[1] << " " << local.start[2] << "; "
              << local.end[0] << " " << local.end[1] << " " << local.end[2] << std::endl;

    stk::parallel_machine_finalize();
    return 0;
}
