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
    //! Type used for sizing of arrays
    using size_t = decltype(sizeof 1);

    //! Index types which can have signed values
    using idx_t = int;

    //! Default dimension of the structured grid
    static constexpr int ndim = 3;
};

/** A structured cartesian aligned box within the index space
 *
 *  Represents an arbitrary box in the index space. The extents are given by
 *  `start` and `end` arrays, The `end` array has the same meaning is one past
 *  the last valid index, such that `size = end - start`. If the box contains
 *  ghosts, then the number of ghosts are represented by the `nghost` array
 *  along each dimension. The bounding box includes the ghosts along each
 *  dimension.
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

/** Structured Cartesian mesh
 *
 *  A structured mesh that has the global information represented by the
 *  `global` StructBox and the local partition in `local` instance object. In
 *  addition, it also carries information of the partitioning when used on
 *  multiple MPI ranks.
 *
 *  The global grid is always defined by the size of the mesh (IMAX, JMAX,
 *  KMAX), and with starting indices always set to 0 and has 0 ghost layers. The
 *  local partition may contain ghost layers.
 */
class StructGrid
{
public:
    StructGrid() = default;

    ~StructGrid() = default;

    /**
     *  @param nghost number of ghost cell layers
     */
    StructGrid(int nghost)
        : nghost_{nghost, nghost, nghost}
    {}

    /** Initialize the global grid extents
     *
     *  Takes three arguments, the number of cells along each dimension that is
     *  used to set the size of the global mesh.
     */
    void set_global_grid(SGTraits::idx_t, SGTraits::idx_t, SGTraits::idx_t);

    /** Decompose the mesh into MPI partitions
     *
     *  The mesh is subdivided into nblocks = (px * py * pz) and each block is
     *  assigned to a unique MPI rank. If `nblocks` is not equal to the total
     *  MPI_Comm size, then it this method will throw an error. This method will
     *  also update the `local` instance with the parition info for the mesh
     *  block residing in the current MPI rank.
     *
     *  If ghosting is desired in the local block, then set_num_ghost must be
     *  called prior to initializing the partitions.
     */
    void set_partitions(int px = 1, int py = 1, int pz = 1);

    /** Set the number of ghost cell layers for sub-blocks in each MPI rank
     */
    void set_num_ghost(int nghost)
    {
        for (int i=0; i < SGTraits::ndim; ++i)
            nghost_[i] = nghost;
    }

    //! Return number of ghost cell layers for this grid
    const int* num_ghosts() const { return nghost_; }

    //! The StructBox instance that has information about the global mesh
    const StructBox& global() const { return global_; }

    //! The StructBox instance that has information about the local block on
    //! this MPI rank
    const StructBox& local() const { return local_; }

    //! Return the number of blocks along each dimension
    const int* partitions() const { return partitions_; }

    //! Return the index of the local block on this MPI rank along each
    //! dimension
    const int* partindex() const { return partindex_; }

private:
    //! Global extents of the grid across all MPI ranks
    StructBox global_;

    //! Local grid info
    StructBox local_;

    //! Partitions along each dimension
    int partitions_[SGTraits::ndim];

    //! Index of this local block
    int partindex_[SGTraits::ndim];

    //! Ghost layers along each direction
    int nghost_[SGTraits::ndim] = {0, 0, 0};

    //! Flag indicating whether the global grid has been initialized
    bool gridInitialized_{false};

    //! Flag indicating whether the local partition has been initialized
    bool partInitialized_{false};
};

}  // nalu
}  // sierra


#endif /* STRUCTGRID_H */
