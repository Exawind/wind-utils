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

#ifndef STRUCTGRIDUTILS_H
#define STRUCTGRIDUTILS_H

#include <stdexcept>
#include <iostream>

#include "struct_grid/StructGrid.h"
#include "struct_grid/StructGridIx.h"
#include "struct_grid/FVStencil.h"

namespace sierra {
namespace nalu {

namespace sgix {

/** Return a box that contains the ghost layer in the given direction
 *
 *  This function computes the box for only the ghost layer in the given
 *  direction that overlaps with the "real box" for the other directions.
 *
 *  @param box The input box for which ghost layer is computed
 *  @param dir The direction of the desired ghost layer (FVStencil neighborhoods)
 *
 *  @return Box containing the indices of the ghost layer
 */
inline StructBox
ghost_layer(const StructBox& box, const fvm::FVStencil dir)
{
    StructBox gbox = real_box(box);

    switch (dir) {
    case fvm::WEST:
        gbox.start[0] = box.start[0];
        gbox.end[0] = box.start[0] + box.nghost[0];
        gbox.size[0] = gbox.end[0] - gbox.start[0];
        break;

    case fvm::EAST:
        gbox.end[0] = box.end[0];
        gbox.start[0] = box.end[0] - box.nghost[0];
        gbox.size[0] = gbox.end[0] - gbox.start[0];
        break;

    case fvm::SOUTH:
        gbox.start[1] = box.start[1];
        gbox.end[1] = box.start[1] + box.nghost[1];
        gbox.size[1] = gbox.end[1] - gbox.start[1];
        break;

    case fvm::NORTH:
        gbox.end[1] = box.end[1];
        gbox.start[1] = box.end[1] - box.nghost[1];
        gbox.size[1] = gbox.end[1] - gbox.start[1];
        break;

    case fvm::BOTTOM:
        gbox.start[2] = box.start[2];
        gbox.end[2] = box.start[2] + box.nghost[2];
        gbox.size[2] = gbox.end[2] - gbox.start[2];
        break;

    case fvm::TOP:
        gbox.end[2] = box.end[2];
        gbox.start[2] = box.end[2] - box.nghost[2];
        gbox.size[2] = gbox.end[2] - gbox.start[2];
        break;

    default:
        throw std::runtime_error("Invalid ghost direction specified");
    }

    return gbox;
}

} // namespace sgix
}  // nalu
}  // sierra


#endif /* STRUCTGRIDUTILS_H */
