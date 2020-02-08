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

#include <numeric>

#include "HexBlockBase.h"
#include "core/PerfUtils.h"
#include "core/ParallelInfo.h"
#include "struct_grid/StructGridIx.h"

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
    const auto& pinfo = get_mpi();
    const std::string timerName("HexBlockBase::initialize");
    auto timeMon = get_stopwatch(timerName);

    pinfo.info() << "HexBlockBase: Registering parts to meta data" << std::endl;

    VectorFieldType& coords = meta_.declare_field<VectorFieldType>(
        stk::topology::NODE_RANK, "coordinates");
    stk::mesh::put_field_on_mesh(coords, meta_.universal_part(), meta_.spatial_dimension(), nullptr);

    // Block mesh part
    {
        stk::mesh::Part& part = meta_.declare_part(
            blockName_, stk::topology::ELEM_RANK);
        stk::io::put_io_part_attribute(part);
        stk::mesh::set_topology(part, stk::topology::HEX_8);

        pinfo.info() << "\tMesh block: " << blockName_ << std::endl;

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

    stk::parallel_machine_barrier(get_mpi().comm());
}

void HexBlockBase::generate_elements()
{
    using EntID = stk::mesh::EntityId;
    const auto& pinfo = get_mpi();

    stk::mesh::Part* part = (meta_.get_part(blockName_));
    stk::mesh::PartVector meshParts{part};

    const auto& global = elemGrid_.global();
    const int* gsize = global.size;
    EntID numGlobalElems = 1;
    EntID numGlobalNodes = 1;
    EntID numNodes = 1;
    EntID numElems = 1;
    for (int i=0; i<3; i++) {
        numNodes *= static_cast<EntID>(meshDims_[i] + 1);
        numElems *= static_cast<EntID>(meshDims_[i]);
        numGlobalNodes *= static_cast<EntID>(gsize[i] + 1);
        numGlobalElems *= static_cast<EntID>(gsize[i]);
    }

    pinfo.info() << "Num. nodes = " << numGlobalNodes << "; Num elements = "
                 << numGlobalElems << std::endl;

    std::vector<stk::mesh::EntityId> nodeIDs(numNodes), elemIDs(numElems);
    stk::mesh::EntityVector nodes, elems;
    nodes.reserve(numNodes);
    elems.reserve(numElems);

    EntID mx = meshDims_[0];
    EntID my = meshDims_[1];
    EntID mz = meshDims_[2];
    EntID nx = meshDims_[0] + 1;
    EntID ny = meshDims_[1] + 1;

    std::vector<stk::mesh::EntityId> nids(8);
    bulk_.modification_begin("Adding mesh nodes");
    {
        {
            pinfo.info() << "\tGenerating nodes..." << std::flush;
            const std::string timerName("HexBlockBase::create_nodes");
            auto timeMon = get_stopwatch(timerName);

            const EntID start =
                nodeIDStart_ + (
                    static_cast<EntID>(nodeBlock_.start[2]) *
                    static_cast<EntID>(nodeBlock_.size[0]) *
                    static_cast<EntID>(nodeBlock_.size[1]));
            std::iota(nodeIDs.begin(), nodeIDs.end(), start);
            ThrowRequire(nodeIDs[0] == start);
            ThrowRequire(nodeIDs[0] < numGlobalNodes);
            ThrowRequire(nodeIDs[numNodes - 1] <= numGlobalNodes);
            bulk_.declare_entities(
                stk::topology::NODE_RANK, nodeIDs, meshParts, nodes);

            if (pinfo.size() > 1) {
                if (pinfo.rank() > 0)
                    add_node_sharing(nodes, 0, pinfo.rank() - 1);
                if (pinfo.rank() < (pinfo.size() - 1))
                    add_node_sharing(
                        nodes, (nodeBlock_.size[2] - 1), (pinfo.rank() + 1));
            }
            pinfo.info() << "done" << std::endl;
            stk::parallel_machine_barrier(pinfo.comm());
        }

        {
            const std::string timerName("HexBlockBase::create_elements");
            auto timeMon = get_stopwatch(timerName);
            pinfo.info() << "\tGenerating elements..." << std::flush;
            const auto& local = elemGrid_.local();
            const EntID start =
                elemIDStart_ + (
                    static_cast<EntID>(local.start[2]) *
                    static_cast<EntID>(local.size[0]) *
                    static_cast<EntID>(local.size[1]));
            std::iota(elemIDs.begin(), elemIDs.end(), start);
            ThrowRequire(elemIDs[0] == start);
            ThrowRequire(elemIDs[0] < numGlobalElems);
            ThrowRequire(elemIDs[numElems - 1] <= numGlobalElems);
            bulk_.declare_entities(
                stk::topology::ELEM_RANK, elemIDs, meshParts, elems);
            pinfo.info() << "done" << std::endl;
        }

        {
            const std::string timerName("HexBlockBase::create_connectivity");
            auto timeMon = get_stopwatch(timerName);
            auto perm = stk::mesh::Permutation::INVALID_PERMUTATION;
            stk::mesh::OrdinalVector scratch1, scratch2, scratch3;

            pinfo.info() << "\tCreating element connectivity... " << std::flush;
            EntID idx;
            for (EntID k=0; k < mz; k++) {
                EntID ik = k * (nx * ny);
                EntID ikp1 = (k+1) * (nx * ny);
                for (EntID j=0; j < my; j++) {
                    EntID ij = j * nx;
                    EntID ijp1 = (j+1) * nx;
                    for (EntID i=0; i < mx; i++) {
                        idx = k*(mx*my) + j*mx + i;
                        nids[0] = ik + ij + i;
                        nids[1] = ik + ij + i + 1;
                        nids[2] = ik + ijp1 + i + 1;
                        nids[3] = ik + ijp1 + i;

                        nids[4] = ikp1 + ij + i;
                        nids[5] = ikp1 + ij + i + 1;
                        nids[6] = ikp1 + ijp1 + i + 1;
                        nids[7] = ikp1 + ijp1 + i;

                        for (EntID ni=0; ni < 8; ++ni)
                            bulk_.declare_relation(
                                elems[idx], nodes[nids[ni]], ni, perm, scratch1,
                                scratch2, scratch3);
                    }
                }
            }
            pinfo.info() << "done" << std::endl;
            stk::parallel_machine_barrier(pinfo.comm());
        }

        {
            const std::string timerName("HexBlockBase::create_sidesets");
            auto timeMon = get_stopwatch(timerName);
            {
                generate_x_boundary(elems, XMIN);
                generate_x_boundary(elems, XMAX);
                generate_y_boundary(elems, YMIN);
                generate_y_boundary(elems, YMAX);
                generate_z_boundary(elems, ZMIN);
                generate_z_boundary(elems, ZMAX);
                stk::parallel_machine_barrier(pinfo.comm());
            }
        }
    }
    stk::parallel_machine_barrier(pinfo.comm());
    pinfo.info() << "\tFinalizing bulk data modifications ... " << std::flush;
    bulk_.modification_end();
    pinfo.info() << "done" << std::endl;

    elemIDs.clear();
    nodes.clear();
    elems.clear();

    pinfo.info() << "\tGenerating coordinates..." << std::endl;
    generate_coordinates(nodeIDs);
    stk::parallel_machine_barrier(pinfo.comm());
}

void HexBlockBase::generate_x_boundary(
    const std::vector<stk::mesh::Entity>& elemVec,
    const SideIDType id)
{
    using EntID = stk::mesh::EntityId;
    EntID mx = meshDims_[0];
    EntID my = meshDims_[1];
    EntID mz = meshDims_[2];

    unsigned sideOrd = 0;
    std::string ssname;
    int ix = 0;
    get_sideset_info(id, ix, ssname, sideOrd);
    get_mpi().info() << "\tGenerating X Sideset: " << ssname << std::endl;
    stk::mesh::Part* part = meta_.get_part(ssname);
    stk::mesh::PartVector partVec{part};

    for (EntID k=0; k<mz; k++) {
        for (EntID j=0; j<my; j++) {
            const EntID idx = k * (mx * my) + j * mx + ix;
            bulk_.declare_element_side(elemVec[idx], sideOrd, partVec);
        }
    }
}

void HexBlockBase::generate_y_boundary(
    const std::vector<stk::mesh::Entity>& elemVec,
    const SideIDType id)
{
    using EntID = stk::mesh::EntityId;
    EntID mx = meshDims_[0];
    EntID my = meshDims_[1];
    EntID mz = meshDims_[2];

    unsigned sideOrd = 0;
    std::string ssname;
    int iy = 0;
    get_sideset_info(id, iy, ssname, sideOrd);
    get_mpi().info() << "\tGenerating Y Sideset: " << ssname << std::endl;
    stk::mesh::Part* part = meta_.get_part(ssname);
    stk::mesh::PartVector partVec{part};

    for (EntID k=0; k<mz; k++) {
        for (EntID i=0; i<mx; i++) {
            const EntID idx = k * (mx * my) + iy * mx + i;
            bulk_.declare_element_side(elemVec[idx], sideOrd, partVec);
        }
    }
}

void HexBlockBase::generate_z_boundary(
    const std::vector<stk::mesh::Entity>& elemVec,
    const SideIDType id)
{
    using EntID = stk::mesh::EntityId;
    const auto& pinfo = get_mpi();
    if (!(((id == ZMIN) && (pinfo.rank() == 0)) ||
          ((id == ZMAX) && (pinfo.rank() == (pinfo.size() - 1)))))
        return;

    EntID mx = meshDims_[0];
    EntID my = meshDims_[1];

    unsigned sideOrd = 0;
    std::string ssname;
    int iz = 0;
    get_sideset_info(id, iz, ssname, sideOrd);
    pinfo.info() << "\tGenerating Z Sideset: " << ssname << std::endl;
    stk::mesh::Part* part = meta_.get_part(ssname);
    stk::mesh::PartVector partVec{part};

    for (EntID j=0; j<my; j++) {
        for (EntID i=0; i<mx; i++) {
            const EntID idx = iz * (mx * my) + j * mx + i;
            bulk_.declare_element_side(elemVec[idx], sideOrd, partVec);
        }
    }
}

void
HexBlockBase::add_node_sharing(
    const std::vector<stk::mesh::Entity>& nodes,
    const unsigned zindex,
    const int otherProc)
{
    using EntID = stk::mesh::EntityId;
    const EntID nx = meshDims_[0] + 1;
    const EntID ny = meshDims_[1] + 1;
    const EntID offset = zindex * nx * ny;

    EntID idx = 0;
    for (EntID j=0; j < ny; ++j)
        for (EntID i = 0; i < nx; ++i) {
            bulk_.add_node_sharing(nodes[offset + idx], otherProc);
            idx++;
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
