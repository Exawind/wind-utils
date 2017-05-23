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
 * With the above input definition, it will use the bounding box of the
 * `fluid_part` to determine the bounding box of the plane to be generated. This
 * will provide coordinate axis aligned sapling planes in x and y directions.
 * Alternately, the user can specify `boundary_type` to be `quad_vertices` and
 * provide the vertices of the quadrilateral that is used to generate the
 * sampling plane as shown below:
 *
 * ```
 * generate_planes:
 *   boundary_type: quad_vertices
 *   fluid_part: Unspecified-2-hex
 *
 *   heights: [ 50.0, 70.0, 90.0 ]
 *   part_name_format: "zplane_%06.1f"
 *
 *   nx: 25  # Number of divisions along (1-2) and (4-3) vertices
 *   ny: 25  # Number of divisions along (1-4) and (2-3) vertices
 *   vertices:
 *     - [250.0, 0.0]
 *     - [500.0, -250.0]
 *     - [750.0, 0.0]
 *     - [500.0, 250.0]
 * ```
 *
 * `part_name_format` is a printf-like format specification that takes one
 * argument - the height as a floating point number. The user can use this to
 * tailor how the nodesets or the shell parts are named in the output Exodus
 * file.
 *
 */
class SamplingPlanes: public PreProcessingTask
{
public:
    /** Sampling Plane boundary type
     */
    enum PlaneBoundaryType {
        BOUND_BOX = 0, ///< Use bounding box of the fluid mesh defined
        QUAD_VERTICES  ///< Use user-defined vertex list for plane boundary
    };

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

    std::vector<std::vector<double>> vertices_;

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

    //! User defined selection of plane boundary type
    PlaneBoundaryType bdyType_{BOUND_BOX};
};

}  // nalu
}  // sierra

#endif /* SAMPLINGPLANES_H */
