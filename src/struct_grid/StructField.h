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

#ifndef STRUCTFIELD_H
#define STRUCTFIELD_H

#include <vector>

#include "struct_grid/StructGrid.h"
#include "struct_grid/StructGridIx.h"

namespace sierra {
namespace nalu {

template<typename T, typename Indexer>
class StructBoxField
{
public:
    using idx_t = SGTraits::idx_t;

    StructBoxField(const StructBox& box)
        : bx_(box),
          ix_(bx_),
          field_(sgix::num_cells(box))
    {}

    inline T& operator()(idx_t i, idx_t j, idx_t k)
    {
        return field_[ix_(i, j, k)];
    }

    inline const T& operator()(idx_t i, idx_t j, idx_t k) const
    {
        return field_[ix_(i, j, k)];
    }

    inline T& operator[](int idx)
    {
        return field_[idx];
    }

    inline const T& operator[](int idx) const
    {
        return field_[idx];
    }

private:
    const StructBox bx_;

    const Indexer ix_;

    std::vector<T> field_;
};

template<typename T, typename Layout = sgix::LeftLayout>
using BoxField = StructBoxField<T, sgix::Indexer<Layout>>;

template <typename T, typename Layout = sgix::LeftLayout>
using PeriodicBoxField = StructBoxField<T, sgix::PeriodicIndexer<Layout>>;

}  // nalu
}  // sierra


#endif /* STRUCTFIELD_H */
