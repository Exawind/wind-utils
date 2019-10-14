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

#ifndef HEXBLOCKBASE_H
#define HEXBLOCKBASE_H

#include "core/CFDMesh.h"
#include "core/YamlUtils.h"
#include "core/ClassRegistry.h"
#include "struct_grid/StructGrid.h"

#include <vector>

namespace sierra {
namespace nalu {

/** Base class representation of a structured hex mesh
 */
class HexBlockBase
{
public:
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

    HexBlockBase(CFDMesh&);

    virtual ~HexBlockBase();

    /** Registers the element block and the sidesets to the STK MetaData instance
     */
    virtual void initialize();

    /** Creates the nodes and elements within the mesh block, processes
     * sidesets, and initializes the coordinates of the mesh structure.
     */
    virtual void run();

    DECLARE_INHERITANCE_REGISTRY
    (
        HexBlockBase,
        (
            CFDMesh& mesh,
            const YAML::Node& node
        ),
        (mesh, node)
    );

    /** Runtime creation of mesh generator instance
     */
    static HexBlockBase* create(
        CFDMesh&,
        const YAML::Node&,
        std::string);

private:
    HexBlockBase() = delete;
    HexBlockBase(const HexBlockBase&) = delete;

protected:
    //! Load common input parameters
    virtual void load(const YAML::Node&);

    //! Generate elements
    virtual void generate_elements();

    //! Generate the xmin and xmax sidesets
    virtual void generate_x_boundary(const std::vector<stk::mesh::Entity> &,
                                     const SideIDType);

    //! Generate the ymin and ymax sidesets
    virtual void generate_y_boundary(const std::vector<stk::mesh::Entity> &,
                                     const SideIDType);

    //! Generate the zmin and zmax sidesets
    virtual void generate_z_boundary(const std::vector<stk::mesh::Entity> &,
                                     const SideIDType);
    //! Generate coordinates
    virtual void generate_coordinates(const std::vector<stk::mesh::EntityId>&) = 0;

    /** Sideset information helper
     * \param[in] id Sideset min/max id
     * \param[inout] index Sideset coordinate index
     * \param[inout] name Sideset name
     * \param[inout] ord Sideset number
     */
    virtual void get_sideset_info(const SideIDType, int&, std::string&, unsigned&);

    virtual void add_node_sharing(
        const std::vector<stk::mesh::Entity>&, const unsigned, const int);

    //! STK Metadata object
    stk::mesh::MetaData& meta_;

    //! STK Bulkdata object
    stk::mesh::BulkData& bulk_;

    //! Block information in parallel runs
    StructGrid elemGrid_;

    StructBox nodeBlock_;

    stk::mesh::EntityId nodeIDStart_{1};

    stk::mesh::EntityId elemIDStart_{1};

    //! Mesh dimensions in each direction
    std::vector<int> meshDims_;

    //! Name of the fluid domain block
    std::string blockName_{"fluid"};

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

};

}  // nalu
}  // sierra



#endif /* HEXBLOCKBASE_H */
