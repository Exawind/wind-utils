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

#ifndef STRUCTGRID_H
#define STRUCTGRID_H

namespace sierra {
namespace nalu {

struct SGTraits
{
    using size_t = decltype(sizeof 1);
    using idx_t = int;
    static constexpr int ndim = 3;
};

/** Representation of a structured box
 */
struct StructBox
{
    //! Dimensions of the box
    static constexpr int ndim = SGTraits::ndim;

    //! Grid sizes in each direction
    SGTraits::idx_t size[ndim] = {0, 0, 0};

    //! Starting index for the box
    SGTraits::idx_t start[ndim] = {0, 0, 0};

    /** The ending index for each dimension in the box
     *
     *  This index follows the one past the last convention, i.e.,
     *  size[i] = end[i] - start[i]
     */
    SGTraits::idx_t end[ndim] = {-1, -1, -1};

    //! Ghost cell layer
    SGTraits::idx_t nghost[ndim] = {0, 0, 0};
};

class StructGrid
{
public:
    StructGrid() = default;

    ~StructGrid() = default;

    StructGrid(int nghost)
        : nghost_{nghost, nghost, nghost}
    {}

    void set_global_grid(SGTraits::idx_t, SGTraits::idx_t, SGTraits::idx_t);

    void set_partitions(int px = 1, int py = 1, int pz = 1);

    void set_num_ghost(int nghost)
    {
        for (int i=0; i < SGTraits::ndim; ++i)
            nghost_[i] = nghost;
    }

    const int* num_ghosts() const { return nghost_; }

    const StructBox& global() const { return global_; }

    const StructBox& local() const { return local_; }

private:
    //! Global extents of the grid across all MPI ranks
    StructBox global_;

    //! Local grid info
    StructBox local_;

    int partitions_[SGTraits::ndim];

    int partindex_[SGTraits::ndim];

    int nghost_[SGTraits::ndim] = {0, 0, 0};

    bool gridInitialized_{false};

    bool partInitialized_{false};
};

}  // nalu
}  // sierra


#endif /* STRUCTGRID_H */
