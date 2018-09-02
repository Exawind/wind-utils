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

#include "BdyIOPlanes.h"
#include "core/PerfUtils.h"
#include "stk_mesh/base/FEMHelpers.hpp"
#include "stk_mesh/base/Bucket.hpp"

namespace sierra {
namespace nalu {

REGISTER_DERIVED_CLASS(PreProcessingTask, BdyIOPlanes, "create_bdy_io_mesh");

BdyIOPlanes::BdyIOPlanes(
    CFDMesh& mesh,
    const YAML::Node& node
) : PreProcessingTask(mesh),
    mesh_(mesh),
    iomesh_(mesh.comm(), mesh_.meta().spatial_dimension()),
    bdyNames_(0)
{
    if (mesh_.meta().spatial_dimension() != 3)
        throw std::runtime_error("Boundary IO plane generation only available for 3-D meshes");
    load(node);
}

BdyIOPlanes::~BdyIOPlanes()
{}

void BdyIOPlanes::load(const YAML::Node& node)
{
    // IO realm mesh
    output_db_ = node["output_db"].as<std::string>();

    // Parse boundary names that needs to be extracted
    const auto& bNames = node["boundary_parts"];
    if (bNames.Type() == YAML::NodeType::Scalar) {
        bdyNames_.push_back(bNames.as<std::string>());
    } else {
        bdyNames_ = bNames.as<std::vector<std::string>>();
    }
}

void BdyIOPlanes::initialize()
{
    const std::string timerName = "BdyIOPlanes::initialize";
    auto timeMon = get_stopwatch(timerName);
    auto& fmeta = mesh_.meta();
    auto& iometa = iomesh_.meta();

    // We need coordinates
    VectorFieldType& coords = iometa.declare_field<VectorFieldType>(
        stk::topology::NODE_RANK, "coordinates");

    for (auto bdyName: bdyNames_) {
        stk::mesh::Part* part = fmeta.get_part(bdyName);
        if (nullptr == part)
            throw std::runtime_error(
                "create_bdy_io_mesh: Invalid boundary part specified = " +
                bdyName);
        stk::mesh::Part& iopart = iometa.declare_part(
            bdyName, stk::topology::ELEM_RANK);
        stk::mesh::set_topology(iopart, stk::topology::SHELL_QUAD_4);
        stk::io::put_io_part_attribute(iopart);
        stk::mesh::put_field_on_mesh(coords, iopart, iometa.spatial_dimension(), nullptr);
    }
}

void BdyIOPlanes::create_boundary(const std::string bdyName)
{
    const std::string timerName = "BdyIOPlanes::create_boundary";
    auto timeMon = get_stopwatch(timerName);
    auto& fmeta = mesh_.meta();
    auto& fbulk = mesh_.bulk();
    auto& iometa = iomesh_.meta();
    auto& iobulk = iomesh_.bulk();
    auto iproc = fbulk.parallel_rank();
    bool doPrint = (fbulk.parallel_rank() == 0);

    if (doPrint)
        std::cout << "\t- " << bdyName << "... ";

    stk::mesh::Part* fpart = fmeta.get_part(bdyName);
    stk::mesh::Part* iopart = iometa.get_part(bdyName);

    stk::mesh::Selector fbdySel = stk::mesh::Selector(*fpart);
    const auto& bkts = fbulk.get_buckets(
        stk::topology::NODE_RANK, fbdySel);

    // Determine the nodes that need to be created on this proc
    size_t num_nodes = 0;
    for (auto b: bkts) num_nodes += b->size();

    std::vector<stk::mesh::EntityId> newIDs(num_nodes);
    std::unordered_map<stk::mesh::EntityId, stk::mesh::EntityId> nodeMap;

    iobulk.modification_begin();
    {
        {
            const std::string timerName = "BdyIOPlanes::create_nodes";
            auto timeMon = get_stopwatch(timerName);
            iobulk.generate_new_ids(stk::topology::NODE_RANK, num_nodes, newIDs);

            size_t nodeIdx=0;
            for (auto b: bkts) {
                for (size_t in=0; in<b->size(); in++) {
                    auto fnode = (*b)[in];
                    auto fnodeid = fbulk.identifier(fnode);
                    auto ionodeid = newIDs[nodeIdx++];
                    auto node = iobulk.declare_entity(stk::topology::NODE_RANK, ionodeid, *iopart);
                    nodeMap[fnodeid] = ionodeid;

                    // TODO: `in_shared` is marked for deprecation. Determine from
                    // STK team what the proposed alternative is.
                    if (fbulk.in_shared(fbulk.entity_key(fnode)))
                        iobulk.add_node_sharing(node, iproc);
                }
            }
        }

        {
            const std::string timerName = "BdyIOPlanes::create_elements";
            auto timeMon = get_stopwatch(timerName);
            const auto& facebkts = fbulk.get_buckets(
                fmeta.side_rank(), (fbdySel & fmeta.locally_owned_part()));

            size_t num_elems = 0;
            for (auto b: facebkts) num_elems += b->size();
            std::vector<stk::mesh::EntityId> elemIDs(num_elems);

            iobulk.generate_new_ids(stk::topology::ELEM_RANK, num_elems, elemIDs);

            size_t eidx=0;
            stk::mesh::EntityIdVector nids(4);
            for (auto b: facebkts) {
                assert(b->topology() == stk::topology::QUAD_4);
                for (size_t k=0; k<b->size(); k++) {
                    const stk::mesh::Entity* fnrels = b->begin_nodes(k);

                    for (int in=0; in < 4; in++)
                        nids[in] = nodeMap[fbulk.identifier(fnrels[in])];

                    stk::mesh::declare_element(iobulk, *iopart, elemIDs[eidx++], nids);
                }
            }
        }
    }
    iobulk.modification_end();

    // Copy coordinates
    VectorFieldType* fcoords = fmeta.get_field<VectorFieldType>(
        stk::topology::NODE_RANK, "coordinates");
    VectorFieldType* iocoords = iometa.get_field<VectorFieldType>(
        stk::topology::NODE_RANK, "coordinates");

    {
        const std::string timerName = "BdyIOPlanes::update_coordinates";
        auto timeMon = get_stopwatch(timerName);
        size_t idx=0;
        for (auto b: bkts) {
            for (size_t i=0; i < b->size(); i++) {
                auto fnode = (*b)[i];
                auto ionode = iobulk.get_entity(stk::topology::NODE_RANK, newIDs[idx++]);

                double* fpt = stk::mesh::field_data(*fcoords, fnode);
                double* iopt = stk::mesh::field_data(*iocoords, ionode);

                for (int j=0; j<3; j++) {
                    iopt[j] = fpt[j];
                }
            }
        }
    }

    if (doPrint)
        std::cout << "done" << std::endl;
}


void BdyIOPlanes::run()
{
    const std::string timerName = "BdyIOPlanes::run";
    auto timeMon = get_stopwatch(timerName);
    bool doPrint = (mesh_.bulk().parallel_rank() == 0);
    iomesh_.meta().commit();

    if (doPrint)
        std::cout << "Extracting boundary planes: " << std::endl;
    for (auto name: bdyNames_)
        create_boundary(name);

    if (doPrint)
        std::cout << "Writing IO mesh: " << output_db_ << std::endl;
    auto& stkio = iomesh_.stkio();
    stkio.set_bulk_data(iomesh_.bulk());
    iomesh_.write_database(output_db_);
}

}  // nalu
}  // sierra
