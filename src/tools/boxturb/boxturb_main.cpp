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
#include <fstream>

#include "tools/boxturb/BoxTurb.h"
#include "tools/boxturb/BoxTurbIO.h"
#include "core/YamlUtils.h"
#include "core/ParallelInfo.h"

#include "boost/program_options.hpp"
#include "stk_util/parallel/Parallel.hpp"

int main(int argc, char**argv)
{
    stk::ParallelMachine comm = stk::parallel_machine_init(&argc, &argv);
    const auto& pinfo = sierra::nalu::get_mpi();

    std::string filename;
    boost::program_options::options_description desc(
        "Nalu turbsim convertor utility. Valid options are");
    desc.add_options()("help,h", "Show this help message")(
        "input-file,i",
        boost::program_options::value<std::string>(&filename)->default_value(
            "boxturb.yaml"),
        "Input file with preprocessor options");

    boost::program_options::variables_map vmap;
    boost::program_options::store(
        boost::program_options::parse_command_line(argc, argv, desc), vmap);

    if (vmap.count("help")) {
        pinfo.info() << desc << std::endl;
        return 0;
    }

    filename = vmap["input-file"].as<std::string>();
    std::ifstream fin(filename.c_str());
    if (!fin.good()) {
        pinfo.info() << "Cannot find input file: " << filename << std::endl;
        return 1;
    }

    pinfo.info() << "\nNalu Turbulent File Processing Utility\n"
                 << "Input file: " << filename << std::endl;

    YAML::Node inpfile(YAML::LoadFile(filename));
    const auto& node = inpfile["boxturb"];

    sierra::nalu::BoxTurb turb;
    turb.load(node);
    turb.run(node);

    stk::parallel_machine_finalize();
    return 0;
}
