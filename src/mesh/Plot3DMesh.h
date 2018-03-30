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

#ifndef PLOT3DMESH_H
#define PLOT3DMESH_H

#include "HexBlockBase.h"

namespace sierra {
namespace nalu {

class Plot3DMesh : public HexBlockBase
{
public:
    Plot3DMesh(
        CFDMesh&,
        const YAML::Node&);

    virtual ~Plot3DMesh();

private:
    Plot3DMesh() = delete;
    Plot3DMesh(const Plot3DMesh&) = delete;

    //! Process the YAML input data and initialize class members
    void load(const YAML::Node&);

    //! Generate coordinates
    void generate_coordinates(const std::vector<stk::mesh::EntityId>&);

    //! Parse Plot3D headers
    void parse_p3d_headers();

    std::string p3dFile_;

    //! Total length of the Plot3D header entry in bytes
    int skipBytes_{20};
};

}  // nalu
}  // sierra


#endif /* PLOT3DMESH_H */
