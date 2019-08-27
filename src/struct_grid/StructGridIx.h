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

#ifndef STRUCTGRIDIX_H
#define STRUCTGRIDIX_H

#include <type_traits>

#include "struct_grid/StructGrid.h"

namespace sierra {
namespace nalu {
namespace sgix {

struct LeftLayout {};
struct RightLayout {};

inline SGTraits::size_t
num_cells(const StructBox& box)
{
    SGTraits::size_t counter = 1;
    for (int i = 0; i < SGTraits::ndim; ++i)
        counter *= box.size[i];

    return counter;
}

inline StructBox
grow_box(
    const StructBox& box,
    SGTraits::idx_t ngx,
    SGTraits::idx_t ngy,
    SGTraits::idx_t ngz)
{
    StructBox box1;
    SGTraits::idx_t ghosts[SGTraits::ndim] = {ngx, ngy, ngz};

    for (int i = 0; i < SGTraits::ndim; ++i) {
        box1.nghost[i] = box.nghost[i] + ghosts[i];
        box1.start[i]  = box.start[i] - ghosts[i];
        box1.end[i]    = box.end[i] + ghosts[i];
        box1.size[i]   = box1.end[i] - box1.start[i];
    }

    return box1;
}

inline StructBox
grow_box(const StructBox& box, SGTraits::idx_t nghost = 1)
{
    return grow_box(box, nghost, nghost, nghost);
}

inline StructBox
real_box(const StructBox& box)
{
    StructBox box1;

    for (int i = 0; i < SGTraits::ndim; ++i) {
        box1.start[i] = box.start[i] + box.nghost[i];
        box1.end[i] = box.end[i] - box.nghost[i];
        box1.size[i] = box1.end[i] - box1.start[i];
    }

    return box1;
}

template<typename Layout = LeftLayout>
class Indexer
{
public:
    using idx_t = SGTraits::idx_t;

    Indexer(const StructBox& box) : bx(box) {}

    template<typename U = Layout>
    inline SGTraits::size_t operator()(
        idx_t i, idx_t j, idx_t k,
        typename std::enable_if<std::is_same<U, LeftLayout>::value, int>::
            type = 0) const
    {
        return (((i - bx.start[0]) * bx.size[1] * bx.size[2]) +
                ((j - bx.start[1]) * bx.size[2]) + (k - bx.start[2]));
    }

    template<typename U = Layout>
    inline SGTraits::size_t operator()(
        idx_t i, idx_t j, idx_t k,
        typename std::enable_if<std::is_same<U, RightLayout>::value, int>::
            type = 0) const
    {
        return (((k - bx.start[2]) * bx.size[0] * bx.size[1]) +
                ((j - bx.start[1]) * bx.size[0]) + (i - bx.start[0]));
    }

private:
    const StructBox& bx;
};

template<typename Layout = LeftLayout>
class PeriodicIndexer
{
public:
    using idx_t = SGTraits::idx_t;

    PeriodicIndexer(const StructBox& box) : bx(box) {}

    template<typename U = Layout>
    inline SGTraits::size_t operator()(
        idx_t i, idx_t j, idx_t k,
        typename std::enable_if<std::is_same<U, LeftLayout>::value, int>::
            type = 0) const
    {
        idx_t ijk[SGTraits::ndim] = {i, j, k};

        for (int d=0; d < SGTraits::ndim; ++d) {
            if (ijk[d] < bx.start[d])
                ijk[d] = bx.end[d] - (bx.start[d] - ijk[d]);
            else if (ijk[d] >= bx.end[d])
                ijk[d] = bx.start[d] + (ijk[d] - bx.end[d]);
        }

        return (((ijk[0] - bx.start[0]) * bx.size[1] * bx.size[2]) +
                ((ijk[1] - bx.start[1]) * bx.size[2]) +
                 (ijk[2] - bx.start[2]));
    }

    template<typename U = Layout>
    inline SGTraits::size_t operator()(
        idx_t i, idx_t j, idx_t k,
        typename std::enable_if<std::is_same<U, RightLayout>::value, int>::
            type = 0) const
    {
        idx_t ijk[SGTraits::ndim] = {i, j, k};

        for (int d = 0; d < SGTraits::ndim; ++d) {
            if (ijk[d] < bx.start[d])
                ijk[d] = bx.end[d] - (bx.start[d] - ijk[d]);
            else if (ijk[d] >= bx.end[d])
                ijk[d] = bx.start[d] + (ijk[d] - bx.end[d]);
        }

        return (((ijk[2] - bx.start[2]) * bx.size[0] * bx.size[1]) +
                ((ijk[1] - bx.start[1]) * bx.size[0]) + (ijk[0] - bx.start[0]));
    }

private:
    const StructBox& bx;
};

template <typename LambdaFunction>
void
ijk_loop(const StructBox& box, LambdaFunction func)
{
    using idx_t = SGTraits::idx_t;
    for (idx_t i = box.start[0]; i < box.end[0]; ++i)
        for (idx_t j = box.start[1]; j < box.end[1]; ++j)
            for (idx_t k = box.start[2]; k < box.end[2]; ++k)
                func(i, j, k);
}

template <typename LambdaFunction>
void
kji_loop(const StructBox& box, LambdaFunction func)
{
    using idx_t = SGTraits::idx_t;
    for (idx_t k = box.start[2]; k < box.end[2]; ++k)
        for (idx_t j = box.start[1]; j < box.end[1]; ++j)
            for (idx_t i = box.start[0]; i < box.end[0]; ++i)
                func(i, j, k);
}

} // namespace sgix
} // namespace nalu
} // namespace sierra

#endif /* STRUCTGRIDIX_H */
