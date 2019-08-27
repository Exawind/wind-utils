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

#ifndef NEIGHBORMAP_H
#define NEIGHBORMAP_H

#include <stdexcept>

#include "struct_grid/StructGrid.h"
#include "struct_grid/StructGridIx.h"
#include "struct_grid/FVStencil.h"

namespace sierra {
namespace nalu {

template<typename Layout = sgix::RightLayout>
class NeighborMap
{
public:
    NeighborMap(const StructGrid& grid);

    ~NeighborMap() = default;

    inline int operator()(fvm::FVStencil) const;

private:
    const StructGrid& grid_;

    StructBox map_;

    const sgix::PeriodicIndexer<Layout> indexer_;
};

template <typename Layout>
NeighborMap<Layout>::NeighborMap(const StructGrid& grid)
    : grid_(grid), indexer_(map_)
{
    const auto* partitions = grid_.partitions();
    for (int d = 0; d < SGTraits::ndim; ++d) {
        map_.end[d] = partitions[d];
        map_.size[d] = partitions[d];
    }
}

template<typename Layout>
int NeighborMap<Layout>::operator()(fvm::FVStencil dir) const
{
    int proc = -1;
    const auto* pid = grid_.partindex();

    switch (dir) {
    case fvm::WEST:
        proc = indexer_(pid[0] - 1, pid[1], pid[2]);
        break;

    case fvm::EAST:
        proc = indexer_(pid[0] + 1, pid[1], pid[2]);
        break;

    case fvm::SOUTH:
        proc = indexer_(pid[0], pid[1] - 1, pid[2]);
        break;

    case fvm::NORTH:
        proc = indexer_(pid[0], pid[1] + 1, pid[2]);
        break;

    case fvm::BOTTOM:
        proc = indexer_(pid[0], pid[1], pid[2] - 1);
        break;

    case fvm::TOP:
        proc = indexer_(pid[0], pid[1], pid[2] + 1);
        break;

    default:
        throw std::runtime_error("Invalid neighbor request");
    }

    return proc;
}

}  // nalu
}  // sierra


#endif /* NEIGHBORMAP_H */
