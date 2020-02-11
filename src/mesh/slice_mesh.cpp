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
 *  Slice mesh generation utility
 */

#include "mesh/Slice.h"
#include "core/YamlUtils.h"
#include "core/PerfUtils.h"
#include "core/ParallelInfo.h"

#include "stk_util/parallel/Parallel.hpp"
#include "stk_util/environment/OptionsSpecification.hpp"
#include "stk_util/environment/ParseCommandLineArgs.hpp"
#include "stk_io/WriteMesh.hpp"
#include "Kokkos_Core.hpp"

#include <iostream>
#include <fstream>
#include <memory>
#include <string>

int main(int argc, char** argv)
{
    stk::ParallelMachine comm = stk::parallel_machine_init(&argc, &argv);
    const auto& pinfo = sierra::nalu::get_mpi();

    Kokkos::initialize(argc, argv);
    {
        std::string filename;
        stk::OptionsSpecification desc(
            "Nalu turbsim convertor utility. Valid options are");
        desc.add_options()
            ("help,h", "Show this help message")
            ("input-file,i", "Input file for slice mesh",
             stk::TargetPointer<std::string>(&filename),
             stk::DefaultValue<std::string>("slice_mesh.yaml"));

        stk::ParsedOptions vmap;
        stk::parse_command_line_args(
            argc, const_cast<const char**>(argv), desc, vmap);

        if (vmap.count("help")) {
            pinfo.info() << desc << std::endl;
            return 0;
        }

        std::ifstream fin(filename.c_str());
        if (!fin.good()) {
            pinfo.info() << "Cannot find input file: " << filename << std::endl;
            return 1;
        }

        pinfo.info() << "\nSlice Mesh Generation Utility\n"
                     << "Input file: " << filename << std::endl;

        YAML::Node inpfile(YAML::LoadFile(filename));
        const auto& node = inpfile["slice_mesh"];
        const std::string output_db = node["output_db"].as<std::string>();

        const auto& node_slices = node["slices"];
        const int numSliceDefs = node_slices.size();
        std::vector<std::unique_ptr<sierra::nalu::Slice>> slices(numSliceDefs);

        sierra::nalu::CFDMesh mesh(comm, 3);
        pinfo.info() << "Loading slice inputs... " << std::endl;
        for (int i=0; i < numSliceDefs; ++i) {
            slices[i].reset(new sierra::nalu::Slice(mesh));
            slices[i]->load(node_slices[i]);
        }
        pinfo.info() << "Initializing slices... " << std::endl;
        for (int i=0; i < numSliceDefs; ++i) {
            slices[i]->initialize();
        }
        mesh.meta().commit();

        for (int i=0; i < numSliceDefs; ++i) {
            slices[i]->run();
        }
        std::cout << "Writing mesh to file: " << output_db << std::endl;
        auto& stkio = mesh.stkio();
        stkio.set_bulk_data(mesh.bulk());
        mesh.write_database(output_db);

        sierra::nalu::summarize_memory_usage(comm, pinfo.info());
    }
    Kokkos::finalize_all();

    stk::parallel_machine_finalize();
    return 0;
}
