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

/** \file
 *
 * Nalu Preprocessing Utility
 *
 * Usage:
 *    `mpiexec -np 1 nalu_preprocess -i nalu_preprocess.yaml`
 *
 */

#include "preprocessing/PreProcessDriver.h"
#include "core/KokkosWrappers.h"

#include "stk_util/parallel/Parallel.hpp"

#include "stk_util/environment/OptionsSpecification.hpp"
#include "stk_util/environment/ParseCommandLineArgs.hpp"

#include <iostream>
#include <memory>
#include <fstream>

int main(int argc, char** argv)
{
    stk::ParallelMachine comm = stk::parallel_machine_init(&argc, &argv);
    Kokkos::initialize(argc, argv);

    std::string inpfile;
    stk::OptionsSpecification desc(
        "Nalu preprocessor utility. Valid options are");
    desc.add_options()
        ("help,h", "Show this help message")
        ("input-file,i",
         "Input file with preprocessor options",
         stk::TargetPointer<std::string>(&inpfile),
         stk::DefaultValue<std::string>("nalu_preprocess.yaml"));

    stk::ParsedOptions vmap;
    stk::parse_command_line_args(
        argc, const_cast<const char**>(argv), desc, vmap);

    if (vmap.count("help")) {
        if (!stk::parallel_machine_rank(comm))
            std::cout << desc << std::endl;
        return 0;
    }

    std::ifstream fin(inpfile.c_str());
    if (!fin.good()) {
        if (!stk::parallel_machine_rank(comm)) {
            std::cout << "Cannot find input file: " << inpfile << std::endl;
        }
        return 1;
    }

    if (stk::parallel_machine_rank(comm) == 0) {
        std::cout << "\nNalu Preprocessing Utility" << "\n"
                  << "Input file: " << inpfile << std::endl;
    }
    Kokkos::DefaultExecutionSpace::print_configuration(std::cout);
    {
        sierra::nalu::PreProcessDriver preprocess(comm, inpfile);
        preprocess.run();
    }

    Kokkos::finalize_all();

    stk::parallel_machine_finalize();
    return 0;
}
