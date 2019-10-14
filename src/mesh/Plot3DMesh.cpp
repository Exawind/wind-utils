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

#include "Plot3DMesh.h"
#include "core/PerfUtils.h"
#include "core/ParallelInfo.h"

#include "stk_mesh/base/TopologyDimensions.hpp"
#include "stk_mesh/base/FEMHelpers.hpp"
#include "stk_mesh/base/Field.hpp"

#include <fstream>

namespace sierra {
namespace nalu {

REGISTER_DERIVED_CLASS(HexBlockBase, Plot3DMesh, "convert_plot3d");

Plot3DMesh::Plot3DMesh(
    CFDMesh& mesh,
    const YAML::Node& node
) : HexBlockBase(mesh)
{
    if (get_mpi().size() > 1)
        throw std::runtime_error("Plot3DMesh does not support MPI runs");

    load(node);
}

Plot3DMesh::~Plot3DMesh()
{}

void Plot3DMesh::load(const YAML::Node& node)
{
    HexBlockBase::load(node);

    p3dFile_ = node["plot3d_file"].as<std::string>();

    parse_p3d_headers();
}

void Plot3DMesh::parse_p3d_headers()
{
    const bool doPrint = (bulk_.parallel_rank() == 0);
    meshDims_.resize(meta_.spatial_dimension());
    const int ndim = meta_.spatial_dimension();
    int nblocks;
    int intval;

    std::ifstream p3d(p3dFile_, std::ios::in | std::ios::binary);
    if (!p3d.is_open())
        throw std::runtime_error("Plot3DMesh:: Error opening file " + p3dFile_);

    // Check if the number of blocks is provided
    p3d.read(reinterpret_cast<char *>(&intval), sizeof(int));
    if (intval == sizeof(int)) {
      p3d.read(reinterpret_cast<char *>(&nblocks), sizeof(int));
      // We need this to be one block for now
      if (nblocks != 1)
        throw std::runtime_error("Plot3dMesh:: Cannot handle multiple blocks");
      p3d.read(reinterpret_cast<char *>(&intval), sizeof(int));
      p3d.read(reinterpret_cast<char *>(&intval), sizeof(int));
      // Bytes to skip for header data
      skipBytes_ = 64;
    } else {
        skipBytes_ = 20;
    }

    int expected_size = ndim * sizeof(int);
    if (intval != expected_size)
      throw std::runtime_error("Plot3DMesh:: Invalid file metadata");
    p3d.read(reinterpret_cast<char *>(meshDims_.data()), intval);

    if (doPrint)
      std::cout << "Plot3D mesh dimensions: [" << meshDims_[0] << " x "
                << meshDims_[1] << " x " << meshDims_[2] << "]" << std::endl;

    // We keep track of elements and not nodes in meshDims_
    for (int d=0; d < ndim; d++)
        meshDims_[d]--;

    elemGrid_.set_global_grid(meshDims_[0], meshDims_[1], meshDims_[2]);
    elemGrid_.set_partitions(1, 1, 1);
    p3d.close();
}

void Plot3DMesh::generate_coordinates(const std::vector<stk::mesh::EntityId>& nodeVec)
{
    const std::string timerName("Plot3DMesh::generate_coordinates");
    auto timeMon = get_stopwatch(timerName);
    size_t num_nodes = nodeVec.size();
    std::vector<double> buffer(num_nodes);
    const int ndim = meta_.spatial_dimension();
    const int nbytes = num_nodes * sizeof(double);
    const int nx = meshDims_[0] + 1;
    const int ny = meshDims_[1] + 1;
    const int nz = meshDims_[2] + 1;
    int intval;

    VectorFieldType* coords = meta_.get_field<VectorFieldType>(
        stk::topology::NODE_RANK, "coordinates");

    std::ifstream p3d(p3dFile_, std::ios::in | std::ios::binary);
    if (!p3d.is_open())
        throw std::runtime_error("Plot3DMesh:: Error opening file " + p3dFile_);

    // Skip over the header
    p3d.seekg(skipBytes_, std::ios::beg);
    for (int d=0; d < ndim; d++) {
        p3d.read(reinterpret_cast<char *>(&intval), sizeof(int));
        if (intval != nbytes)
            throw std::runtime_error(
                "Plot3DMesh:: Inconsistent data sizes encountered while reading coordinate "
                + std::to_string(d+1));

        p3d.read(reinterpret_cast<char*>(buffer.data()), nbytes);
        p3d.read(reinterpret_cast<char *>(&intval), sizeof(int));

        size_t idx = 0;
        for (int k=0; k < nz; k++)
            for (int j=0; j < ny; j++)
                for (int i=0; i < nx; i++) {
                    auto node = bulk_.get_entity(
                        stk::topology::NODE_RANK, nodeVec[idx]);
                    double* pt = stk::mesh::field_data(*coords, node);
                    pt[d] = buffer[idx++];
                }
    }
}

}  // nalu
}  // sierra
