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

#include "mesh/HexBlockMesh.h"
#include "core/YamlUtils.h"

#include "stk_util/parallel/Parallel.hpp"
#include "stk_io/WriteMesh.hpp"

#include "boost/program_options.hpp"

#include <iostream>
#include <fstream>
#include <memory>

int main(int argc, char** argv)
{
    stk::ParallelMachine comm = stk::parallel_machine_init(&argc, &argv);

    std::string inpfile;
    boost::program_options::options_description desc(
        "Nalu ABL mesh generation utility. Valid options are");
    desc.add_options()
        ("help,h", "Show this help message")
        ("input-file,i",
         boost::program_options::value<std::string>(&inpfile)->default_value(
             "nalu_abl_mesh.yaml"),
         "Input file with preprocessor options");

    boost::program_options::variables_map vmap;
    boost::program_options::store(
        boost::program_options::parse_command_line(argc, argv, desc), vmap);

    if (vmap.count("help")) {
        if (!stk::parallel_machine_rank(comm))
            std::cout << desc << std::endl;
        return 0;
    }

    inpfile = vmap["input-file"].as<std::string>();
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

    YAML::Node doc(YAML::LoadFile(inpfile));
    const YAML::Node node = doc["nalu_abl_mesh"];
    const std::string output_db = node["output_db"].as<std::string>();

    std::unique_ptr<sierra::nalu::CFDMesh> mesh(new sierra::nalu::CFDMesh(comm, 3));

    sierra::nalu::HexBlockMesh blockMesh(*mesh, node);

    blockMesh.initialize();
    mesh->meta().commit();

    blockMesh.run();

    std::cout << "Writing mesh to file: " << output_db << std::endl;
    bool set_64bit = true;
    sierra::nalu::wind_utils::get_optional(node, "ioss_8bit_ints", set_64bit);

    if (set_64bit) mesh->set_64bit_flags();
    auto& stkio = mesh->stkio();
    stkio.set_bulk_data(mesh->bulk());
    mesh->write_database(output_db);

    stk::parallel_machine_finalize();
    return 0;
}
