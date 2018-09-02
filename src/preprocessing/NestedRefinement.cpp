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

#include "NestedRefinement.h"
#include "core/YamlUtils.h"
#include "core/ClassRegistry.h"
#include "core/PerfUtils.h"

namespace {
using sierra::nalu::NestedRefinement;

inline void cross_prod(
    const NestedRefinement::Vec3D& v1,
    const NestedRefinement::Vec3D& v2,
    NestedRefinement::Vec3D& v3)
{
    v3[0] = v1[1] * v2[2] - v1[2] * v2[1];
    v3[1] = v1[2] * v2[0] - v1[0] * v2[2];
    v3[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

inline double dot_prod(
    const NestedRefinement::Vec3D& v1,
    const NestedRefinement::Vec3D& v2)
{
    return (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]);
}

inline void compute_midpoint(
    const stk::mesh::BulkData& bulk,
    const sierra::nalu::VectorFieldType& coords,
    const stk::mesh::Entity& elem,
    NestedRefinement::Vec3D& midPt)
{
    for (int d=0; d < 3; d++)
        midPt[d] = 0.0;

    auto* nodes = bulk.begin_nodes(elem);
    auto num_nodes = bulk.num_nodes(elem);

    for (size_t in=0; in < num_nodes; in++) {
        double* crd = stk::mesh::field_data(coords, nodes[in]);
        for (int d=0; d < 3; d++)
            midPt[d] += crd[d];
    }

    for (int d=0; d < 3; d++)
        midPt[d] /= num_nodes;
}
}

namespace sierra {
namespace nalu {

REGISTER_DERIVED_CLASS(PreProcessingTask, NestedRefinement, "mesh_local_refinement");

NestedRefinement::NestedRefinement(
    CFDMesh& mesh,
    const YAML::Node& node
) : PreProcessingTask(mesh)
{
    load(node);
}

void
NestedRefinement::load(const YAML::Node& node)
{
    fluidPartNames_ = node["fluid_parts"].as<std::vector<std::string>>();
    auto tLocs = node["turbine_locations"].as<std::vector<std::vector<double>>>();
    numTurbines_ = tLocs.size();
    turbineLocs_.resize(numTurbines_);
    for (size_t i=0; i < numTurbines_; i++)
        for (size_t j=0; j < 3; j++)
            turbineLocs_[i][j] = tLocs[i][j];

    {
        auto& tmp = node["turbine_diameters"];
        if (tmp.Type() == YAML::NodeType::Scalar) {
            turbineDia_.resize(numTurbines_);
            const double tdia = tmp.as<double>();
            for (size_t t=0; t < numTurbines_; t++)
                turbineDia_[t] = tdia;
        } else {
            turbineDia_ = node["turbine_diameters"].as<std::vector<double>>();
            if (turbineDia_.size() != numTurbines_)
                throw std::runtime_error(
                    "NestedRefinement:: Turbine diameters don't match number of turbines");
        }
    }

    {
        auto& tmp = node["turbine_heights"];
        if (tmp.Type() == YAML::NodeType::Scalar) {
            turbineHt_.resize(numTurbines_);
            const double tht = tmp.as<double>();
            for (size_t t=0; t < numTurbines_; t++)
                turbineHt_[t] = tht;
        } else {
            turbineHt_ = node["turbine_heights"].as<std::vector<double>>();
            if (turbineHt_.size() != numTurbines_)
                throw std::runtime_error(
                    "NestedRefinement:: Turbine heights don't match number of turbines");
        }
    }

    refineLevels_ = node["refinement_levels"].as<std::vector<std::vector<double>>>();
    numLevels_ = refineLevels_.size();
    for (size_t i=0; i < numLevels_; i++)
        if (refineLevels_[i].size() != 4)
            throw std::runtime_error("NestedRefinement:: Incorrect sizes for input data.");

    wind_utils::get_optional(node, "search_tolerance", searchTol_);
    wind_utils::get_optional(node, "refine_field_name", refineFieldName_);
    wind_utils::get_optional(node, "write_percept_files", writePercept_);
    wind_utils::get_optional(node, "percept_file_prefix", perceptFilePrefix_);

    corners_.resize(numTurbines_ * numLevels_);
    boxLengths_.resize(numTurbines_ * numLevels_);
    boxAxes_.resize(numTurbines_);

    // Parse the box orientation inputs
    const YAML::Node& orient = node["orientation"];
    std::string orientType = "wind_direction";
    wind_utils::get_optional(orient, "type", orientType);
    if (orientType == "wind_direction") {
        windAngle_ = orient["wind_direction"].as<double>();
    } else {
        throw std::runtime_error("Unsupported option for orientation type.");
    }

    // Create local refinement box information
    process_inputs();
}

void
NestedRefinement::initialize()
{
    const std::string timerName = "NestedRefinement::initialize";
    auto timeMon = get_stopwatch(timerName);
    auto& meta = mesh_.meta();
    ScalarFieldType& refiner = meta.declare_field<ScalarFieldType>(
        stk::topology::ELEM_RANK, refineFieldName_);
    fluidParts_.resize(fluidPartNames_.size());
    for (size_t i=0; i < fluidPartNames_.size(); i++) {
        stk::mesh::Part* part = meta.get_part(fluidPartNames_[i]);
        if (part == nullptr)
            throw std::runtime_error("NestedRefinement: Missing fluid part in database: " +
                                     fluidPartNames_[i]);
        else {
            fluidParts_[i] = part;
            stk::mesh::put_field_on_mesh(refiner, *part, 1, nullptr);
        }
    }
    mesh_.add_output_field(refineFieldName_);
}

void
NestedRefinement::run()
{
    const std::string timerName = "NestedRefinement::run";
    auto timeMon = get_stopwatch(timerName);
    auto& meta = mesh_.meta();
    auto& bulk = mesh_.bulk();
    const bool doPrint = (bulk.parallel_rank() == 0);
    const stk::mesh::Selector sel = meta.locally_owned_part()
        & stk::mesh::selectUnion(fluidParts_);
    const auto& bkts = bulk.get_buckets(
        stk::topology::ELEM_RANK, sel);
    auto* coords = meta.get_field<VectorFieldType>(
        stk::topology::NODE_RANK, "coordinates");
    auto* refiner = meta.get_field<ScalarFieldType>(
        stk::topology::ELEM_RANK, refineFieldName_);

    if (doPrint)
        std::cout << "Processing percept field: " << refineFieldName_ << std::endl;

    Vec3D midPt;
    for (auto b: bkts) {
        for (size_t ie=0; ie < b->size(); ie++) {
            auto elem = (*b)[ie];
            compute_midpoint(bulk, *coords, elem, midPt);

            double* refval = stk::mesh::field_data(*refiner, elem);
            refval[0] = std::max(refval[0], compute_refine_fraction(midPt));
        }
    }

    if (writePercept_) {
        write_percept_inputs();
    } else if (doPrint) {
        std::cout << "Skipping percept input file write" << std::endl;
    }
    mesh_.set_write_flag();
}

double
NestedRefinement::compute_refine_fraction(Vec3D& point)
{
    // pt1 - Position vector of the point w.r.t. min corner of the box
    // pt2 - Position vector (pt1) transformed to the box aligned reference frame
    // frac - fraction used for refinement of this particular element
    Vec3D pt1, pt2;
    double frac = 0.0;
    for (size_t it=0; it < numTurbines_; it++) {
        const size_t offset = it * numLevels_;
        for (size_t iz=0; iz < numLevels_; iz++) {
            // Position vector w.r.t. origin of this box
            for (int d = 0; d < 3; d++)
                pt1[d] = point[d] - corners_[offset+iz][d];

            // Transform to local vector in the box aligned axes
            for (int d=0; d < 3; d++)
                pt2[d] = dot_prod(pt1, boxAxes_[it][d]);

            // Check if the point is inside this refinement box
            bool needRefine = true;
            for (int d=0; d < 3; d++) {
                if ((-searchTol_ > pt2[d]) ||
                    (pt2[d] > boxLengths_[offset+iz][d])) {
                    needRefine = false;
                    // Break early if we find that this point is outside the box
                    break;
                }
            }

            if (needRefine) {
                // We want the turbine closest to this point to win out in terms
                // of refinement level
                frac = std::max(
                    frac, static_cast<double>(iz+1)/static_cast<double>(numLevels_));
            } else {
                // No need to check the nested box if the point is outside the bigger box
                break;
            }
        }
    }
    return frac;
}

void
NestedRefinement::process_inputs()
{
    // We create a local reference frame for each nested box refinement such
    // that the principal axes are such that x axis aligns with the direction of
    // wind, the y-axis is perpendicular to the direction of wind and z-axis is
    // aligned with the z-axis in the global reference frame. The origin
    // (corner_) is the lower corner of this bounding box and the lengths are
    // the dimensions of the box in each of the principal axes. The boxAxes_ are
    // unit vectors along each box and are used to transform any point to this
    // box-local axis.
    //

    // Populate the transformation matrix to convert a point into the box
    // aligned coordinate frame. It is assumed that each nested box for a
    // particular turbine is aligned with the others and so we maintain only one
    // local transformation per turbine.
    //
    // TODO: Generalize this for arbitrary orientations other than just wind direction.
    auto ang = windAngle_ * std::acos(-1.0) / 180.0;
    auto cang = - std::cos(ang);
    auto sang = - std::sin(ang);
    for (size_t t=0; t < numTurbines_; t++) {
        // X-vector
        boxAxes_[t][0][0] = sang;
        boxAxes_[t][0][1] = cang;
        boxAxes_[t][0][2] = 0.0;
        // Z-vector
        boxAxes_[t][2][0] = 0.0;
        boxAxes_[t][2][1] = 0.0;
        boxAxes_[t][2][2] = 1.0;
        // Y-axes as cross product of x and z
        cross_prod(boxAxes_[t][2], boxAxes_[t][0], boxAxes_[t][1]);
    }

    // Generate the lower corner of the bounding box and the lengths in each
    // direction of the box.
    for (size_t it=0; it < numTurbines_; it++) {
        const size_t offset = it * numLevels_;
        for (size_t iz=0; iz < numLevels_; iz++) {
            for (int i=0; i < 3; i++) {
                corners_[offset+iz][i] = turbineLocs_[it][i];
                // z coordinate is also shifted by the tower height
                if (i == 2)
                    corners_[offset + iz][i] += turbineHt_[it];
                for (int j=0; j < 3; j++) {
                    int ix = (j == 0)? 0 : (j + 1);
                    corners_[offset+iz][i] -= (
                         boxAxes_[it][j][i] * refineLevels_[iz][ix] * turbineDia_[it]);
                }
            }
            boxLengths_[offset + iz][0] =
                turbineDia_[it] *
                    (refineLevels_[iz][0] + refineLevels_[iz][1]) +
                searchTol_;
            boxLengths_[offset + iz][1] =
                2.0 * turbineDia_[it] * refineLevels_[iz][2] + searchTol_;
            boxLengths_[offset + iz][2] =
                2.0 * turbineDia_[it] * refineLevels_[iz][3] + searchTol_;
        }
    }
}

void
NestedRefinement::write_percept_inputs()
{
    auto iproc = mesh_.bulk().parallel_rank();
    if (iproc != 0) return;

    YAML::Node node, marker;

    node["error_indicator_field"] = refineFieldName_;
    node["do_rebalance"] = true;
    node["max_refinement_level"] = 10;
    marker["refine_fraction"] = 1.0;
    marker["unrefine_fraction"] = 0.0;
    marker["type"] = "fraction";

    std::cout << "Writing percept input files... " << std::endl;
    for (size_t i=0; i < numLevels_; i++) {
      marker["refine_fraction"] =
          0.9 * static_cast<double>(i + 1) / static_cast<double>(numLevels_);
      marker["unrefine_fraction"] =
          0.75 * static_cast<double>(i + 1) / static_cast<double>(numLevels_);
      node["marker"] = marker;
      std::string fname = perceptFilePrefix_ + std::to_string(i + 1) + ".yaml";
      std::cout << "\t" << fname << std::endl;
      std::ofstream adapt(fname.c_str());
      adapt << node << std::endl;
      adapt.close();
    }
    std::cout << "Sample percept command line: \n"
              << "mesh_adapt --refine=DEFAULT --input_mesh=mesh0.e "
              << "--output_mesh=mesh1.e --RAR_info=adapt1.yaml" << std::endl;
}

}  // nalu
}  // sierra
