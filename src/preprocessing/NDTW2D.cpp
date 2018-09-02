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

#include "NDTW2D.h"

namespace sierra {
namespace nalu {

REGISTER_DERIVED_CLASS(PreProcessingTask, NDTW2D, "calc_ndtw2d_deprecated");

NDTW2D::NDTW2D(
    CFDMesh& mesh,
    const YAML::Node& node
) : PreProcessingTask(mesh),
    meta_(mesh.meta()),
    bulk_(mesh.bulk()),
    fluid_parts_(0),
    wall_parts_(0),
    wall_dist_name_("NDTW"),
    ndim_(meta_.spatial_dimension())
{
    // This is a temporary utility that is not scalable
    if (bulk_.parallel_size() > 1)
        throw std::runtime_error("NDTW2D is not a parallel utility");

    std::cerr << "!!!WARNING!!! NDTW2D is a deprecated utility." << std::endl;
    load(node);
}

void NDTW2D::load(const YAML::Node& wdist)
{
    auto fluid_partnames = wdist["fluid_parts"].as<std::vector<std::string>>();
    auto wall_partnames = wdist["wall_parts"].as<std::vector<std::string>>();

    if(wdist["wall_dist_name"]) {
        wall_dist_name_ = wdist["wall_dist_name"].as<std::string>();
    }

    fluid_parts_.resize(fluid_partnames.size());
    wall_parts_.resize(wall_partnames.size());

    for(size_t i=0; i<fluid_partnames.size(); i++) {
        stk::mesh::Part* part = meta_.get_part(fluid_partnames[i]);
        if (NULL == part) {
            throw std::runtime_error("Missing fluid part in mesh database: " +
                                     fluid_partnames[i]);
        } else {
            fluid_parts_[i] = part;
        }
    }
    for(size_t i=0; i<wall_partnames.size(); i++) {
        stk::mesh::Part* part = meta_.get_part(wall_partnames[i]);
        if (NULL == part) {
            throw std::runtime_error("Missing wall part in mesh database: " +
                                     wall_partnames[i]);
        } else {
            wall_parts_[i] = part;
        }
    }
}

void NDTW2D::initialize()
{
    VectorFieldType* coords = meta_.get_field<VectorFieldType>(
        stk::topology::NODE_RANK, "coordinates");
    ScalarFieldType& ndtw = meta_.declare_field<ScalarFieldType>(
        stk::topology::NODE_RANK, wall_dist_name_);

    for(auto part: fluid_parts_) {
        stk::mesh::put_field_on_mesh(*coords, *part, ndim_, nullptr);
        stk::mesh::put_field_on_mesh(ndtw, *part, nullptr);
    }

}

void NDTW2D::run()
{
    calc_ndtw();
    // Register this field for output during write
    mesh_.add_output_field(wall_dist_name_);

    mesh_.set_write_flag();
}

void NDTW2D::calc_ndtw()
{
    stk::mesh::Selector fluid_union = stk::mesh::selectUnion(fluid_parts_);
    stk::mesh::Selector wall_union = stk::mesh::selectUnion(wall_parts_);

    const stk::mesh::BucketVector& fluid_bkts = bulk_.get_buckets(
        stk::topology::NODE_RANK, fluid_union);
    const stk::mesh::BucketVector& wall_bkts = bulk_.get_buckets(
        stk::topology::NODE_RANK, wall_union);

    VectorFieldType* coords = meta_.get_field<VectorFieldType>(
        stk::topology::NODE_RANK, "coordinates");
    ScalarFieldType* ndtw = meta_.get_field<ScalarFieldType>(
        stk::topology::NODE_RANK, wall_dist_name_);

    std::cout << "Calculating nearest wall distance... " << std::endl;
    for(size_t ib=0; ib < fluid_bkts.size(); ib++) {
        stk::mesh::Bucket& fbkt = *fluid_bkts[ib];
        double* xyz = stk::mesh::field_data(*coords, fbkt);
        double* wdist = stk::mesh::field_data(*ndtw, fbkt);

        for(size_t in=0; in < fbkt.size(); in++) {
            double min_dist = std::numeric_limits<double>::max();
            for(size_t jb=0; jb < wall_bkts.size(); jb++) {
                stk::mesh::Bucket& wbkt = *wall_bkts[jb];
                double* wxyz = stk::mesh::field_data(*coords, wbkt);

                for(size_t jn=0; jn< wbkt.size(); jn++) {
                    double dist_calc = 0.0;
                    for(int j=0; j<ndim_; j++) {
                        double dst = xyz[in*ndim_+j] - wxyz[jn*ndim_+j];
                        dist_calc += dst * dst;
                    }
                    if (dist_calc < min_dist) min_dist = dist_calc;
                }
            }
            wdist[in] = std::sqrt(min_dist);
        }
    }
}

} // nalu
} // sierra
