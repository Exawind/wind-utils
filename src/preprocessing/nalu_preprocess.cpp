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

#include "stk_util/parallel/Parallel.hpp"

#include "boost/program_options.hpp"

#include <iostream>
#include <memory>
#include <fstream>

int main(int argc, char** argv)
{
    stk::ParallelMachine comm = stk::parallel_machine_init(&argc, &argv);

    std::string inpfile;
    boost::program_options::options_description desc(
        "Nalu preprocessor utility. Valid options are");
    desc.add_options()
        ("help,h", "Show this help message")
        ("input-file,i",
         boost::program_options::value<std::string>(&inpfile)->default_value(
             "nalu_preprocess.yaml"),
         "Input file with preprocessor options");

    boost::program_options::variables_map vmap;
    boost::program_options::store(
        boost::program_options::parse_command_line(argc, argv, desc), vmap);

    if (vmap.count("help")) {
        if (!stk::parallel_machine_rank(comm))
            std::cerr << desc << std::endl;
        return 0;
    }

    inpfile = vmap["input-file"].as<std::string>();
    std::ifstream fin(inpfile.c_str());
    if (!fin.good()) {
        if (!stk::parallel_machine_rank(comm)) {
            std::cerr << "Cannot find input file: " << inpfile << std::endl;
        }
        return 1;
    }

    std::cerr << "\nNalu Preprocessing Utility" << "\n    "
              << "Input file: " << inpfile << std::endl;
    sierra::nalu::PreProcessDriver preprocess(comm, inpfile);
    preprocess.run();

    stk::parallel_machine_finalize();
    return 0;
}
