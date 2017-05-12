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

#ifndef SAMPLINGPLANES_H
#define SAMPLINGPLANES_H

#include "PreProcessingTask.h"

namespace sierra {
namespace nalu {

/** Generate 2-D grids/planes for data sampling
 *
 * Currently only generates horizontal planes at user-defined heights.
 *
 * Requires a section `generate_planes` in the input file within the
 * `nalu_preprocess` section:
 *
 * ```
 * generate_planes:
 *   fluid_part: Unspecified-2-hex
 *
 *   heights: [ 70.0 ]
 *   part_name_format: "zplane_%06.1f"
 *
 *   dx: 12.0
 *   dy: 12.0
 * ```
 *
 * `part_name_format` is a printf-like format specification that takes one
 * argument - the height as a floating point number. The user can use this to
 * tailor how the nodesets or the shell parts are named in the output Exodus
 * file.
 *
 * \todo Handle generation of planes in any direction
 * \todo Enable user option to select node_set or shell topology
 */
class SamplingPlanes: public PreProcessingTask
{
public:
    SamplingPlanes(CFDMesh&, const YAML::Node&);

    virtual ~SamplingPlanes() {}

    virtual void initialize();

    virtual void run();

private:
    SamplingPlanes() = delete;
    SamplingPlanes(const SamplingPlanes&) = delete;

    void load(const YAML::Node&);

    //! Use fluid Realm mesh to estimate the x-y bounding box for the sampling
    //! planes.
    void calc_bounding_box();

    //! Generate entities and update coordinates for a given sampling plane
    void generate_zplane(const double);

    //! STK Metadata object
    stk::mesh::MetaData& meta_;

    //! STK Bulkdata object
    stk::mesh::BulkData& bulk_;

    //! Heights where the averaging planes are generated
    std::vector<double> heights_;

    //! Bounding box of the original mesh
    std::array<std::array<double,3>,2> bBox_;

    //! Format specification for the part name
    std::string name_format_;

    //! Fluid realm parts (to determine mesh bounding box)
    std::vector<std::string> fluidPartNames_;

    //! Parts of the fluid mesh (to determine mesh bounding box)
    stk::mesh::PartVector fluidParts_;

    //! Spatial resolution in x and y directions
    double dx_;

    //! Spatial resolution in x and y directions
    double dy_;

    //! Number of nodes in x and y directions
    size_t nx_, ny_;

    //! Number of elements in x and y directions
    size_t mx_, my_;

    //! Dimensionality of the mesh
    int ndim_;
};

}  // nalu
}  // sierra

#endif /* SAMPLINGPLANES_H */
