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

#include "HexBlockMesh.h"
#include "core/YamlUtils.h"

#include "stk_mesh/base/TopologyDimensions.hpp"
#include "stk_mesh/base/FEMHelpers.hpp"
#include "stk_mesh/base/Field.hpp"

namespace sierra{
namespace nalu {

bool check_input_coords(
    std::vector<std::vector<double>> crdVec,
    size_t ncoords,
    size_t ndim=3)
{
    bool coordsOK = true;

    if (crdVec.size() != ncoords) {
        coordsOK = false;
        return coordsOK;
    }

    for (auto crds: crdVec)
        if (crds.size() != ndim) {
            coordsOK = false;
            break;
        }
    return coordsOK;
}

void bbox_to_vertices(
    std::vector<std::vector<double>>& vertices,
    std::vector<std::vector<double>>& bbox)
{
    // S-W corner (0 vertex)
    vertices[0][0] = bbox[0][0];
    vertices[0][1] = bbox[0][1];
    vertices[0][2] = bbox[0][2];
    // S-E corner (1 vertex)
    vertices[1][0] = bbox[1][0];
    vertices[1][1] = bbox[0][1];
    vertices[1][2] = bbox[0][2];
    // N-E corner (2 vertex)
    vertices[2][0] = bbox[1][0];
    vertices[2][1] = bbox[1][1];
    vertices[2][2] = bbox[0][2];
    // N-W corner (3 vertex)
    vertices[3][0] = bbox[0][0];
    vertices[3][1] = bbox[1][1];
    vertices[3][2] = bbox[0][2];

    // S-W corner (4 vertex)
    vertices[4][0] = bbox[0][0];
    vertices[4][1] = bbox[0][1];
    vertices[4][2] = bbox[1][2];
    // S-E corner (5 vertex)
    vertices[5][0] = bbox[1][0];
    vertices[5][1] = bbox[0][1];
    vertices[5][2] = bbox[1][2];
    // N-E corner (6 vertex)
    vertices[6][0] = bbox[1][0];
    vertices[6][1] = bbox[1][1];
    vertices[6][2] = bbox[1][2];
    // N-W corner (7 vertex)
    vertices[7][0] = bbox[0][0];
    vertices[7][1] = bbox[1][1];
    vertices[7][2] = bbox[1][2];
}

HexBlockMesh::HexBlockMesh(
    CFDMesh& mesh,
    const YAML::Node& node
) : meta_(mesh.meta()),
    bulk_(mesh.bulk()),
    vertices_(8, std::vector<double>(3, 0.0))
{
    if (meta_.spatial_dimension() != 3)
        throw std::runtime_error("HexBlockMesh is only valid for 3-D meshes");

    load(node);
}

HexBlockMesh::~HexBlockMesh()
{}

void HexBlockMesh::load(const YAML::Node& node)
{
    using namespace sierra::nalu::wind_utils;

    if (node["spec_type"]) {
        std::string spec_type = node["spec_type"].as<std::string>();

        if (spec_type == "bounding_box")
            vertexDef_ = BOUND_BOX;
        else if (spec_type == "vertices")
            vertexDef_ = VERTICES;
        else
            throw std::runtime_error("HexBlockMesh: Invalid option for spec_type = "
                                     + spec_type);
    }

    if (vertexDef_ == BOUND_BOX) {
        auto bbox = node["vertices"].as<std::vector<std::vector<double>>>();
        auto coordsOK = check_input_coords(bbox, 2, 3);
        if (!coordsOK)
            throw std::runtime_error("HexBlockMesh: Inconsistent coordinates provided");

        bbox_to_vertices(vertices_, bbox);

    } else {
        vertices_ = node["vertices"].as<std::vector<std::vector<double>>>();

        auto coordsOK = check_input_coords(vertices_, 8, 3);
        if (!coordsOK)
            throw std::runtime_error("HexBlockMesh: Inconsistent coordinates provided");
    }

    meshDims_ = node["mesh_dimensions"].as<std::vector<int>>();

    get_optional(node, "fluid_part_name", blockName_);
    get_optional(node, "xmin_boundary_name", ss_xmin_name_);
    get_optional(node, "xmax_boundary_name", ss_xmax_name_);
    get_optional(node, "ymin_boundary_name", ss_ymin_name_);
    get_optional(node, "ymax_boundary_name", ss_ymax_name_);
    get_optional(node, "zmin_boundary_name", ss_zmin_name_);
    get_optional(node, "zmax_boundary_name", ss_zmax_name_);
}

void HexBlockMesh::initialize()
{
    const auto iproc = bulk_.parallel_rank();
    const bool doPrint = (iproc == 0);

    if (doPrint)
        std::cout << "HexBlockMesh: Registering parts to meta data" << std::endl;

    VectorFieldType& coords = meta_.declare_field<VectorFieldType>(
        stk::topology::NODE_RANK, "coordinates");
    stk::mesh::put_field(coords, meta_.universal_part(), meta_.spatial_dimension());

    // Block mesh part
    {
        stk::mesh::Part& part = meta_.declare_part(
            blockName_, stk::topology::ELEM_RANK);
        stk::io::put_io_part_attribute(part);
        stk::mesh::set_topology(part, stk::topology::HEX_8);

        if (doPrint)
            std::cout << "\tMesh block: " << blockName_ << std::endl;

        stk::mesh::put_field(coords, part, meta_.spatial_dimension());
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

void HexBlockMesh::run()
{
    generate_elements();
}

void HexBlockMesh::generate_elements()
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
        if (doPrint)
            std::cout << "\tGenerating node IDs..." << std::endl;

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

        if (doPrint)
            std::cout << "\tGenerating element IDs..." << std::endl;

        bulk_.generate_new_ids(stk::topology::ELEM_RANK, numElems, elemIDs);

        if (doPrint)
            std::cout << "\tCreating elements... " ;
        marker = 1;
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
    bulk_.modification_end();
    if (doPrint)
        std::cout << "\tFinalizing bulk modifications..." << std::endl;

    bulk_.modification_begin("Adding sidesets");
    {
        generate_x_boundary(elemIDs, 0);
        generate_x_boundary(elemIDs, mx-1);
        generate_y_boundary(elemIDs, 0);
        generate_y_boundary(elemIDs, my-1);
        generate_z_boundary(elemIDs, 0);
        generate_z_boundary(elemIDs, mz-1);
    }
    bulk_.modification_end();


    elemIDs.clear();
    elemIDs.resize(1);

    if (doPrint)
        std::cout << "\tGenerating coordinates..." << std::endl;
    generate_coordinates(nodeIDs);
}

void HexBlockMesh::generate_coordinates(const std::vector<stk::mesh::EntityId>& nodeVec)
{
    int nx = meshDims_[0] + 1;
    int ny = meshDims_[1] + 1;
    int nz = meshDims_[2] + 1;

    VectorFieldType* coords = meta_.get_field<VectorFieldType>(
        stk::topology::NODE_RANK, "coordinates");

    // TODO: implement stretching factors
    double dx = 1.0 / static_cast<double>(meshDims_[0]);
    double dy = 1.0 / static_cast<double>(meshDims_[1]);
    double dz = 1.0 / static_cast<double>(meshDims_[2]);

    for (int k=0; k < nz; k++) {
        double rz = k * dz;
        for (int j=0; j < ny; j++) {
            double ry = j * dy;
            for (int i=0; i < nx; i++) {
                double rx = i * dx;
                int idx = k * (nx * ny) + j * nx + i;
                auto node = bulk_.get_entity(stk::topology::NODE_RANK, nodeVec[idx]);
                double* pt = stk::mesh::field_data(*coords, node);

                pt[0] = ((1.0 - rx) * (1.0 - ry) * (1.0 - rz) * vertices_[0][0] +
                         rx * (1.0 - ry) * (1.0 - rz) * vertices_[1][0] +
                         rx * ry * (1.0 - rz) * vertices_[2][0] +
                         (1.0 - rx) * ry * (1.0 - rz) * vertices_[3][0] +
                         (1.0 - rx) * (1.0 - ry) * (rz) * vertices_[4][0] +
                         rx * (1.0 - ry) * (rz) * vertices_[5][0] +
                         rx * ry * (rz) * vertices_[6][0] +
                         (1.0 - rx) * ry * (rz) * vertices_[7][0]);

                pt[1] = ((1.0 - rx) * (1.0 - ry) * (1.0 - rz) * vertices_[0][1] +
                         rx * (1.0 - ry) * (1.0 - rz) * vertices_[1][1] +
                         rx * ry * (1.0 - rz) * vertices_[2][1] +
                         (1.0 - rx) * ry * (1.0 - rz) * vertices_[3][1] +
                         (1.0 - rx) * (1.0 - ry) * (rz) * vertices_[4][1] +
                         rx * (1.0 - ry) * (rz) * vertices_[5][1] +
                         rx * ry * (rz) * vertices_[6][1] +
                         (1.0 - rx) * ry * (rz) * vertices_[7][1]);

                pt[2] = ((1.0 - rx) * (1.0 - ry) * (1.0 - rz) * vertices_[0][2] +
                         rx * (1.0 - ry) * (1.0 - rz) * vertices_[1][2] +
                         rx * ry * (1.0 - rz) * vertices_[2][2] +
                         (1.0 - rx) * ry * (1.0 - rz) * vertices_[3][2] +
                         (1.0 - rx) * (1.0 - ry) * (rz) * vertices_[4][2] +
                         rx * (1.0 - ry) * (rz) * vertices_[5][2] +
                         rx * ry * (rz) * vertices_[6][2] +
                         (1.0 - rx) * ry * (rz) * vertices_[7][2]);
            }
        }
    }
}

void HexBlockMesh::generate_x_boundary(
    const std::vector<stk::mesh::EntityId>& elemVec,
    const int ix)
{
    int mx = meshDims_[0];
    int my = meshDims_[1];
    int mz = meshDims_[2];

    unsigned sideOrd = (ix == 0)? 3 : 1;
    std::string ssname = (ix==0)? ss_xmin_name_ : ss_xmax_name_;
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

void HexBlockMesh::generate_y_boundary(
    const std::vector<stk::mesh::EntityId>& elemVec,
    const int iy)
{
    int mx = meshDims_[0];
    int my = meshDims_[1];
    int mz = meshDims_[2];

    unsigned sideOrd = (iy == 0)? 0 : 2;
    std::string ssname = (iy==0)? ss_ymin_name_ : ss_ymax_name_;
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

void HexBlockMesh::generate_z_boundary(
    const std::vector<stk::mesh::EntityId>& elemVec,
    const int iz)
{
    int mx = meshDims_[0];
    int my = meshDims_[1];

    unsigned sideOrd = (iz == 0)? 4 : 5;
    std::string ssname = (iz==0)? ss_zmin_name_ : ss_zmax_name_;
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

}  // nalu
}  // sierra
