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

#include "HexBlockBase.h"
#include "core/PerfUtils.h"

#include "stk_mesh/base/TopologyDimensions.hpp"
#include "stk_mesh/base/FEMHelpers.hpp"
#include "stk_mesh/base/Field.hpp"

namespace sierra {
namespace nalu {

DEFINE_INHERITANCE_REGISTRY(HexBlockBase);

HexBlockBase*
HexBlockBase::create(
    CFDMesh& mesh,
    const YAML::Node& node,
    std::string lookup)
{
    auto it = HexBlockBaseReg_ConstructorTable_->find(lookup);
    if (it != HexBlockBaseReg_ConstructorTable_->end()) {
        return (it->second)(mesh, node);
    } else {
        std::cout << "ERROR:: Invalid Mesh generation type: " << lookup
                  << std::endl;
        std::cout << "Valid generation types are: " << std::endl;
        for (const auto& t: *HexBlockBaseReg_ConstructorTable_) {
            std::cout << "\t" << t.first << std::endl;
        }
    }
    return nullptr;
}

HexBlockBase::HexBlockBase(
    CFDMesh& mesh
) : meta_(mesh.meta()),
    bulk_(mesh.bulk())
{
    if (meta_.spatial_dimension() != 3)
        throw std::runtime_error("HexBlockBase is only valid for 3-D meshes");
}

HexBlockBase::~HexBlockBase()
{}

void HexBlockBase::load(const YAML::Node& node)
{
    using namespace sierra::nalu::wind_utils;

    get_optional(node, "fluid_part_name", blockName_);
    get_optional(node, "xmin_boundary_name", ss_xmin_name_);
    get_optional(node, "xmax_boundary_name", ss_xmax_name_);
    get_optional(node, "ymin_boundary_name", ss_ymin_name_);
    get_optional(node, "ymax_boundary_name", ss_ymax_name_);
    get_optional(node, "zmin_boundary_name", ss_zmin_name_);
    get_optional(node, "zmax_boundary_name", ss_zmax_name_);
}

void HexBlockBase::initialize()
{
    const auto iproc = bulk_.parallel_rank();
    const bool doPrint = (iproc == 0);
    const std::string timerName("HexBlockBase::initialize");
    auto timeMon = get_stopwatch(timerName);

    if (doPrint)
        std::cout << "HexBlockBase: Registering parts to meta data" << std::endl;

    VectorFieldType& coords = meta_.declare_field<VectorFieldType>(
        stk::topology::NODE_RANK, "coordinates");
    stk::mesh::put_field_on_mesh(coords, meta_.universal_part(), meta_.spatial_dimension(), nullptr);

    // Block mesh part
    {
        stk::mesh::Part& part = meta_.declare_part(
            blockName_, stk::topology::ELEM_RANK);
        stk::io::put_io_part_attribute(part);
        stk::mesh::set_topology(part, stk::topology::HEX_8);

        if (doPrint)
            std::cout << "\tMesh block: " << blockName_ << std::endl;

        stk::mesh::put_field_on_mesh(coords, part, meta_.spatial_dimension(), nullptr);
    }

    // West
    {
        stk::mesh::Part& part = meta_.declare_part(
            ss_xmin_name_, meta_.side_rank());
        stk::io::put_io_part_attribute(part);
    }

    // East
    {
        stk::mesh::Part& part = meta_.declare_part(
            ss_xmax_name_, meta_.side_rank());
        stk::io::put_io_part_attribute(part);
    }

    // South
    {
        stk::mesh::Part& part = meta_.declare_part(
            ss_ymin_name_, meta_.side_rank());
        stk::io::put_io_part_attribute(part);
    }

    // North
    {
        stk::mesh::Part& part = meta_.declare_part(
            ss_ymax_name_, meta_.side_rank());
        stk::io::put_io_part_attribute(part);
    }

    // Terrain
    {
        stk::mesh::Part& part = meta_.declare_part(
            ss_zmin_name_, meta_.side_rank());
        stk::io::put_io_part_attribute(part);
    }

    // Top
    {
        stk::mesh::Part& part = meta_.declare_part(
            ss_zmax_name_, meta_.side_rank());
        stk::io::put_io_part_attribute(part);
    }
}

void HexBlockBase::run()
{
    const std::string timerName("HexBlockBase::run");
    auto timeMon = get_stopwatch(timerName);
    generate_elements();
}

void HexBlockBase::generate_elements()
{
    const auto iproc = bulk_.parallel_rank();
    const bool doPrint = (iproc == 0);

    stk::mesh::Part& part = *(meta_.get_part(blockName_));

    unsigned numNodes = 1;
    unsigned numElems = 1;
    for (int i=0; i<3; i++) {
        numNodes *= (meshDims_[i] + 1);
        numElems *= meshDims_[i];
    }

    if (doPrint) {
        std::cout << "Num. nodes = " << numNodes << "; Num elements = "
                  << numElems << std::endl;
    }

    std::vector<stk::mesh::EntityId> nodeIDs(numNodes), elemIDs(numElems);

    int mx = meshDims_[0];
    int my = meshDims_[1];
    int mz = meshDims_[2];
    int nx = meshDims_[0] + 1;
    int ny = meshDims_[1] + 1;

    stk::mesh::EntityIdVector nids(8);
    bulk_.modification_begin("Adding mesh nodes");
    {
        {
            if (doPrint)
                std::cout << "\tGenerating node IDs..." << std::endl;
            const std::string timerName("HexBlockBase::create_nodes");
            auto timeMon = get_stopwatch(timerName);

            bulk_.generate_new_ids(stk::topology::NODE_RANK, numNodes, nodeIDs);
            if (doPrint)
                std::cout << "\tCreating nodes... " ;
            unsigned marker = 1;
            for (unsigned i=0; i<numNodes; i++) {
                if (doPrint && (marker <= (i * 10 / numNodes))) {
                    std::cout << marker * 10 << "% ";
                    marker++;
                }
                bulk_.declare_entity(
                    stk::topology::NODE_RANK, nodeIDs[i], part);
            }
            std::cout << std::endl;
        }

        {
            const std::string timerName("HexBlockBase::create_elements");
            auto timeMon = get_stopwatch(timerName);
            if (doPrint)
                std::cout << "\tGenerating element IDs..." << std::endl;

            bulk_.generate_new_ids(stk::topology::ELEM_RANK, numElems, elemIDs);

            if (doPrint)
                std::cout << "\tCreating elements... " ;
            unsigned marker = 1;
            int idx;
            for (int k=0; k < mz; k++) {
                int ik = k * (nx * ny);
                int ikp1 = (k+1) * (nx * ny);
                for (int j=0; j < my; j++) {
                    int ij = j * nx;
                    int ijp1 = (j+1) * nx;
                    for (int i=0; i < mx; i++) {
                        idx = k*(mx*my) + j*mx + i;
                        nids[0] = nodeIDs[ik + ij + i];
                        nids[1] = nodeIDs[ik + ij + i + 1];
                        nids[2] = nodeIDs[ik + ijp1 + i + 1];
                        nids[3] = nodeIDs[ik + ijp1 + i];

                        nids[4] = nodeIDs[ikp1 + ij + i];
                        nids[5] = nodeIDs[ikp1 + ij + i + 1];
                        nids[6] = nodeIDs[ikp1 + ijp1 + i + 1];
                        nids[7] = nodeIDs[ikp1 + ijp1 + i];

                        stk::mesh::declare_element(
                            bulk_, part, elemIDs[idx], nids);

                        if (doPrint && (marker <= (idx * 10 / numElems))) {
                            std::cout << marker * 10 << "% ";
                            marker++;
                        }
                    }
                }
            }
            std::cout << std::endl;
        }

    }
    bulk_.modification_end();
    if (doPrint)
        std::cout << "\tFinalizing bulk modifications..." << std::endl;

    {
        const std::string timerName("HexBlockBase::create_sidesets");
        auto timeMon = get_stopwatch(timerName);
        bulk_.modification_begin("Adding sidesets");
        {
            generate_x_boundary(elemIDs, XMIN);
            generate_x_boundary(elemIDs, XMAX);
            generate_y_boundary(elemIDs, YMIN);
            generate_y_boundary(elemIDs, YMAX);
            generate_z_boundary(elemIDs, ZMIN);
            generate_z_boundary(elemIDs, ZMAX);
        }
        bulk_.modification_end();
    }


    elemIDs.clear();
    elemIDs.resize(1);

    if (doPrint)
        std::cout << "\tGenerating coordinates..." << std::endl;
    generate_coordinates(nodeIDs);
}

void HexBlockBase::generate_x_boundary(
    const std::vector<stk::mesh::EntityId>& elemVec,
    const SideIDType id)
{
    int mx = meshDims_[0];
    int my = meshDims_[1];
    int mz = meshDims_[2];

    unsigned sideOrd = 0;
    std::string ssname;
    int ix = 0;
    get_sideset_info(id, ix, ssname, sideOrd);
    std::cout << "\tGenerating X Sideset: " << ssname << std::endl;
    stk::mesh::Part* part = meta_.get_part(ssname);
    stk::mesh::PartVector partVec{part};

    for (int k=0; k<mz; k++) {
        for (int j=0; j<my; j++) {
            const int idx = k * (mx * my) + j * mx + ix;
            auto elem = bulk_.get_entity(stk::topology::ELEM_RANK, elemVec[idx]);
            bulk_.declare_element_side(elem, sideOrd, partVec);
        }
    }
}

void HexBlockBase::generate_y_boundary(
    const std::vector<stk::mesh::EntityId>& elemVec,
    const SideIDType id)
{
    int mx = meshDims_[0];
    int my = meshDims_[1];
    int mz = meshDims_[2];

    unsigned sideOrd = 0;
    std::string ssname;
    int iy = 0;
    get_sideset_info(id, iy, ssname, sideOrd);
    std::cout << "\tGenerating Y Sideset: " << ssname << std::endl;
    stk::mesh::Part* part = meta_.get_part(ssname);
    stk::mesh::PartVector partVec{part};

    for (int k=0; k<mz; k++) {
        for (int i=0; i<mx; i++) {
            const int idx = k * (mx * my) + iy * mx + i;
            auto elem = bulk_.get_entity(stk::topology::ELEM_RANK, elemVec[idx]);
            bulk_.declare_element_side(elem, sideOrd, partVec);
        }
    }
}

void HexBlockBase::generate_z_boundary(
    const std::vector<stk::mesh::EntityId>& elemVec,
    const SideIDType id)
{
    int mx = meshDims_[0];
    int my = meshDims_[1];

    unsigned sideOrd = 0;
    std::string ssname;
    int iz = 0;
    get_sideset_info(id, iz, ssname, sideOrd);
    std::cout << "\tGenerating Z Sideset: " << ssname << std::endl;
    stk::mesh::Part* part = meta_.get_part(ssname);
    stk::mesh::PartVector partVec{part};

    for (int j=0; j<my; j++) {
        for (int i=0; i<mx; i++) {
            const int idx = iz * (mx * my) + j * mx + i;
            auto elem = bulk_.get_entity(stk::topology::ELEM_RANK, elemVec[idx]);
            bulk_.declare_element_side(elem, sideOrd, partVec);
        }
    }
}

void
HexBlockBase::get_sideset_info(
    const SideIDType id, int& index, std::string& name, unsigned& ord)
{
    int mx = meshDims_[0];
    int my = meshDims_[1];
    int mz = meshDims_[2];

    switch (id) {
    case XMIN:
        ord = 3;
        name = ss_xmin_name_;
        index = 0;
        break;
    case XMAX:
        ord = 1;
        name = ss_xmax_name_;
        index = mx - 1;
        break;
    case YMIN:
        ord = 0;
        name = ss_ymin_name_;
        index = 0;
        break;
    case YMAX:
        ord = 2;
        name = ss_ymax_name_;
        index = my - 1;
        break;
    case ZMIN:
        ord = 4;
        name = ss_zmin_name_;
        index = 0;
        break;
    case ZMAX:
        ord = 5;
        name = ss_zmax_name_;
        index = mz - 1;
        break;
    default:
        throw std::runtime_error(
            "Boundary id does not match min/max id value.");
    }
}
}  // nalu
}  // sierra
