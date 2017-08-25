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

#ifndef HEXBLOCKMESH_H
#define HEXBLOCKMESH_H

#include "core/CFDMesh.h"

#include "yaml-cpp/yaml.h"

#include <vector>
#include <array>

namespace sierra {
namespace nalu {

/** Create a structured block mesh with HEX-8 elements
 *
 */
class HexBlockMesh
{
public:
    /** Computational domain definition type
     */
    enum DomainExtentsType {
        BOUND_BOX = 0, ///< Use bounding box to define mesh extents
        VERTICES       ///< Provide vertices for the cuboidal domain
    };

    HexBlockMesh(CFDMesh&, const YAML::Node&);

    virtual ~HexBlockMesh();

    virtual void initialize();

    virtual void run();

private:
    HexBlockMesh() = delete;
    HexBlockMesh(const HexBlockMesh&) = delete;

    //! Process the YAML input data and initialize class members
    void load(const YAML::Node&);

    //! Generate elements
    void generate_elements();

    //! Generate coordinates
    void generate_coordinates(const std::vector<stk::mesh::EntityId>&);

    void generate_x_boundary(const std::vector<stk::mesh::EntityId>&, const int);
    void generate_y_boundary(const std::vector<stk::mesh::EntityId>&, const int);
    void generate_z_boundary(const std::vector<stk::mesh::EntityId>&, const int);

    //! STK Metadata object
    stk::mesh::MetaData& meta_;

    //! STK Bulkdata object
    stk::mesh::BulkData& bulk_;

    //! Corners of the computational domain mesh
    std::vector<std::vector<double>> vertices_;

    //! Mesh dimensions in each direction
    std::array<int, 3> meshDims_;

    //! Stretch factor
    std::array<double, 3> stretchFactor{{1.0, 1.0, 1.0}};

    //! Name of the fluid domain block
    std::string blockName_{"fluid_part"};

    // Names for the boundary sidesets

    std::string ss_xmin_name_{"west"};
    std::string ss_xmax_name_{"east"};
    std::string ss_ymin_name_{"south"};
    std::string ss_ymax_name_{"north"};
    std::string ss_zmin_name_{"terrain"};
    std::string ss_zmax_name_{"top"};

    //! Flag indicating user selection of domain extents
    DomainExtentsType vertexDef_{BOUND_BOX};
};

}  // nalu
}  // sierra

#endif /* HEXBLOCKMESH_H */
