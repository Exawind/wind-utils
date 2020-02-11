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
 *  ABL block mesh generation utility
 */

#include "mesh/HexBlockBase.h"
#include "mesh/HexBlockMesh.h"
#include "core/YamlUtils.h"
#include "core/PerfUtils.h"
#include "core/ParallelInfo.h"

#include "stk_util/parallel/Parallel.hpp"
#include "stk_io/WriteMesh.hpp"

#include "stk_util/environment/OptionsSpecification.hpp"
#include "stk_util/environment/ParseCommandLineArgs.hpp"

#include <iostream>
#include <fstream>
#include <memory>

int main(int argc, char** argv)
{
    stk::ParallelMachine comm = stk::parallel_machine_init(&argc, &argv);

    std::string inpfile;
    stk::OptionsSpecification desc("Nalu ABL mesh generation utility.");
    desc.add_options()
        ("help,h", "Show this help message")
        ("input-file,i",
         "Input file with preprocessor options",
         stk::TargetPointer<std::string>(&inpfile),
         stk::DefaultValue<std::string>("nalu_abl_mesh.yaml"));

    stk::ParsedOptions vmap;
    stk::parse_command_line_args(argc, const_cast<const char**>(argv), desc, vmap);

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
        std::cout << "\nNalu ABL Mesh Generation Utility" << "\n"
                  << "Input file: " << inpfile << std::endl;
    }

    {
        YAML::Node doc(YAML::LoadFile(inpfile));
        const YAML::Node node = doc["nalu_abl_mesh"];
        const std::string output_db = node["output_db"].as<std::string>();

        std::unique_ptr<sierra::nalu::CFDMesh> mesh(
            new sierra::nalu::CFDMesh(comm, 3));

        // sierra::nalu::HexBlockMesh blockMesh(*mesh, node);
        std::string mesh_type = "generate_ablmesh";
        if (node["mesh_type"])
            mesh_type = node["mesh_type"].as<std::string>();
        std::unique_ptr<sierra::nalu::HexBlockBase> blockMesh(
            sierra::nalu::HexBlockBase::create(*mesh, node, mesh_type));

        blockMesh->initialize();
        mesh->meta().commit();

        blockMesh->run();

        const auto& pinfo = sierra::nalu::get_mpi();
        pinfo.info() << "Writing mesh to file: " << output_db << std::endl;
        bool set_64bit = false;
        bool set_auto_join = true;
        sierra::nalu::wind_utils::get_optional(
            node, "ioss_8bit_ints", set_64bit);
        sierra::nalu::wind_utils::get_optional(
            node, "auto_join", set_auto_join);

        if (set_64bit) mesh->set_64bit_flags();
        if (set_auto_join) mesh->set_auto_join();
        auto& stkio = mesh->stkio();
        stkio.set_bulk_data(mesh->bulk());
        mesh->write_database(output_db);
        stk::parallel_machine_barrier(comm);

        sierra::nalu::summarize_memory_usage(comm, std::cout);
        bool print_timing_stats = false;
        sierra::nalu::wind_utils::get_optional(node, "print_timing_stats", print_timing_stats);
        if (print_timing_stats) {
            Teuchos::TimeMonitor::summarize(
                std::cout, false, true, false, Teuchos::Union);
        }
    }

    stk::parallel_machine_finalize();
    return 0;
}
