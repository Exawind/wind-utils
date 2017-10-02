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

#ifndef GEOMETRICSTRETCHING_H
#define GEOMETRICSTRETCHING_H

#include "MeshSpacing.h"

namespace sierra {
namespace nalu {

/** Create a mesh spacing distribution with a constant stretching factor
 *
 *  Requires user to specify a constant stretching factor that is used, along
 *  with the number of elements, to determine the first cell height and the
 *  resulting spacing distribution on a one-dimensional mesh of unit length.
 *  Given a stretching factor \f$s\f$, the first cell height is calculated as
 *
 *  \f[
 *    h_0 = L \left(\frac{s - 1}{s^n - 1}\right)
 *  \f]
 *
 *  By default, the stretching factor is applied in one direction. The user can
 *  set the `bidirectional` flag to true to apply the stretching factors and
 *  spacings at both ends.
 */
class GeometricStretching : public MeshSpacing
{
public:
    GeometricStretching(
        int,
        const YAML::Node&);

    // Initialze a geometrically stretched spacing for a 1-D mesh
    virtual void init_spacings();

private:
    GeometricStretching() = delete;
    GeometricStretching(const GeometricStretching&) = delete;

    //! Process user-provided YAML inputs
    void load(const YAML::Node&);

    //! Helper function to initialize unidirectional spacing distribution
    void unidirectional_stretching();

    //! Helper function to initialize bidirectional spacing distribution
    void bidirectional_stretching();

    //! First cell height (non-dimensional)
    double fch_{0.0};

    //! Stretching factor
    double fac_{1.0};

    //! Flag indicating whether the geometric stretching is applied at both ends.
    bool bidirectional_{false};

    //! Flag indicating whether to dump spacing array to screen
    bool verbose_{false};
};

}  // nalu
}  // sierra


#endif /* GEOMETRICSTRETCHING_H */
