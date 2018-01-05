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
#include "spacing/MeshSpacing.h"

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

    /** Sideset definition type
     */
    enum SideIDType {
        XMIN = 0,
        YMIN,
        ZMIN,
        XMAX,
        YMAX,
        ZMAX
    };

    /**
     * \param mesh A sierra::nalu::CFDMesh instance
     * \param node The YAML::Node containing inputs for this task
     */
    HexBlockMesh(CFDMesh&, const YAML::Node&);

    virtual ~HexBlockMesh();

    /** Registers the element block and the sidesets to the STK MetaData instance
     */
    virtual void initialize();

    /** Creates the nodes and elements within the mesh block, processes
     * sidesets, and initializes the coordinates of the mesh structure.
     */
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

    //! Generate the xmin and xmax sidesets
    void generate_x_boundary(const std::vector<stk::mesh::EntityId>&, const SideIDType);

    //! Generate the ymin and ymax sidesets
    void generate_y_boundary(const std::vector<stk::mesh::EntityId>&, const SideIDType);

    //! Generate the zmin and zmax sidesets
    void generate_z_boundary(const std::vector<stk::mesh::EntityId>&, const SideIDType);

    /** Sideset information helper
     * \param[in] id Sideset min/max id
     * \param[inout] index Sideset coordinate index
     * \param[inout] name Sideset name
     * \param[inout] ord Sideset number
     */
    void get_sideset_info(const SideIDType id, int& index, std::string& name, unsigned& ord);

    //! STK Metadata object
    stk::mesh::MetaData& meta_;

    //! STK Bulkdata object
    stk::mesh::BulkData& bulk_;

    //! Corners of the computational domain mesh
    std::vector<std::vector<double>> vertices_;

    //! Mesh dimensions in each direction
    std::vector<int> meshDims_;

    std::unique_ptr<MeshSpacing> xSpacing_;
    std::unique_ptr<MeshSpacing> ySpacing_;
    std::unique_ptr<MeshSpacing> zSpacing_;

    //! Name of the fluid domain block
    std::string blockName_{"fluid_part"};

    // Names for the boundary sidesets

    //! Name of the XMIN sideset
    std::string ss_xmin_name_{"west"};

    //! Name of the XMAX sideset
    std::string ss_xmax_name_{"east"};

    //! Name of the YMIN sideset
    std::string ss_ymin_name_{"south"};

    //! Name of the YMAX sideset
    std::string ss_ymax_name_{"north"};

    //! Name of the ZMIN sideset
    std::string ss_zmin_name_{"terrain"};

    //! Name of the ZMAX sideset
    std::string ss_zmax_name_{"top"};

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
