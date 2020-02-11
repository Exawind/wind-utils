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

#include "stk_util/environment/OptionsSpecification.hpp"
#include "stk_util/environment/ParseCommandLineArgs.hpp"
#include "stk_util/parallel/Parallel.hpp"

int main(int argc, char**argv)
{
    stk::ParallelMachine comm = stk::parallel_machine_init(&argc, &argv);
    const auto& pinfo = sierra::nalu::get_mpi();

    std::string inpfile;
    stk::OptionsSpecification desc(
        "Nalu postprocessor utility. Valid options are");
    desc.add_options()
        ("help,h", "Show this help message")
        ("input-file,i",
         "Input file with post-processing options",
         stk::TargetPointer<std::string>(&inpfile),
         stk::DefaultValue<std::string>("nalu_postprocess.yaml"));

    stk::ParsedOptions vmap;
    stk::parse_command_line_args(
        argc, const_cast<const char**>(argv), desc, vmap);

    if (vmap.count("help")) {
        pinfo.info() << desc << std::endl;
        return 0;
    }

    std::ifstream fin(inpfile.c_str());
    if (!fin.good()) {
        pinfo.info() << "Cannot find input file: " << inpfile << std::endl;
        return 1;
    }

    pinfo.info() << "\nNalu Turbulent File Processing Utility\n"
                 << "Input file: " << inpfile << std::endl;

    YAML::Node ynode(YAML::LoadFile(inpfile));
    const auto& node = ynode["boxturb"];

    sierra::nalu::BoxTurb turb;
    turb.load(node);
    turb.run(node);

    stk::parallel_machine_finalize();
    return 0;
}
