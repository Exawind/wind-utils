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

#ifndef FVSTENCIL_H
#define FVSTENCIL_H

/** \file FVStencil.h
 *  Second-order finite volume utilities
 */

namespace sierra {
namespace nalu {
namespace fvm {

/** A simple 2nd order finte volume stencil
 */
enum FVStencil
{
    CENTER = 0,
    WEST,
    EAST,
    SOUTH,
    NORTH,
    BOTTOM,
    TOP,
    NUM_STENCIL
};

/** Indices for the neighbors w.r.t. cell center based on FVStencil
 */
static int fv_offsets[NUM_STENCIL][3] = {
    {0, 0, 0},
    {-1, 0, 0},
    {1, 0, 0},
    {0, -1, 0},
    {0, 1, 0},
    {0, 0, -1},
    {0, 0, 1}
};

} // namespace fvm
}  // nalu
}  // sierra


#endif /* FVSTENCIL_H */
