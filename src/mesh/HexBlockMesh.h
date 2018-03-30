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

#include "HexBlockBase.h"
#include "spacing/MeshSpacing.h"

#include "yaml-cpp/yaml.h"

#include <vector>
#include <array>

namespace sierra {
namespace nalu {

/** Create a structured block mesh with HEX-8 elements
 *
 */
class HexBlockMesh : public HexBlockBase
{
public:
    /** Computational domain definition type
     */
    enum DomainExtentsType {
        BOUND_BOX = 0, ///< Use bounding box to define mesh extents
        VERTICES       ///< Provide vertices for the cuboidal domain
    };

    /**
     * \param mesh A sierra::nalu::CFDMesh instance
     * \param node The YAML::Node containing inputs for this task
     */
    HexBlockMesh(CFDMesh&, const YAML::Node&);

    virtual ~HexBlockMesh();

private:
    HexBlockMesh() = delete;
    HexBlockMesh(const HexBlockMesh&) = delete;

    //! Process the YAML input data and initialize class members
    void load(const YAML::Node&);

    //! Generate coordinates
    void generate_coordinates(const std::vector<stk::mesh::EntityId>&);

    //! Corners of the computational domain mesh
    std::vector<std::vector<double>> vertices_;

    std::unique_ptr<MeshSpacing> xSpacing_;
    std::unique_ptr<MeshSpacing> ySpacing_;
    std::unique_ptr<MeshSpacing> zSpacing_;

    //! Spacing type in x-direction
    std::string xspacing_type_{"constant_spacing"};

    //! Spacing type in y-direction
    std::string yspacing_type_{"constant_spacing"};

    //! Spacing type in z-direction
    std::string zspacing_type_{"constant_spacing"};

    //! Flag indicating user selection of domain extents
    DomainExtentsType vertexDef_{BOUND_BOX};
};

}  // nalu
}  // sierra

#endif /* HEXBLOCKMESH_H */
