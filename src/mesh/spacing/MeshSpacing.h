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

#ifndef MESHSPACING_H
#define MESHSPACING_H

#include "core/ClassRegistry.h"

#include "yaml-cpp/yaml.h"

namespace sierra {
namespace nalu {

/** Abstract base class that defines the notion of mesh spacing
 *
 *  This class provides an interface where mesh spacing for a structured mesh
 *  can be represented as a 1-D array of values (`0.0 <= ratio[i] <= 1.0`) in a
 *  particular direction, that represents the location of the i-th node on the
 *  mesh on a unit cube.
 *
 *  \sa sierra::nalu::HexBlockMesh
 */
class MeshSpacing
{
public:
    MeshSpacing(int npts)
        : numPts_(npts),
          ratios_(npts, 0.0)
    {}

    virtual ~MeshSpacing() {}

    //! Initialize spacings based on user inputs
    virtual void init_spacings() = 0;

    //! A 1-D array of fractions that represents the distance from the origin
    //! for a unit cube.
    inline const std::vector<double>& ratios() const
    { return ratios_; }

    DECLARE_INHERITANCE_REGISTRY
    (
        MeshSpacing,
        (
            int npts,
            const YAML::Node& node
        ),
        (npts, node)
    );

    /** Runtime creation of the concrete spacing instance
     */
    static MeshSpacing* create(
        int npts,
        const YAML::Node& node,
        std::string lookup);

private:
    MeshSpacing() = delete;
    MeshSpacing(const MeshSpacing&) = delete;

protected:
    //! Number of points on the mesh in a given direction
    int numPts_;

    //! A 1-D array of fractions that represents the distance from the origin
    //! for a unit cube.
    std::vector<double> ratios_;
};

}  // nalu
}  // sierra


#endif /* MESHSPACING_H */
