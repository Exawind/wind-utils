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

#include "struct_grid/StructGrid.h"
#include "struct_grid/StructGridIx.h"
#include "core/ParallelInfo.h"

namespace sierra {
namespace nalu {

void StructGrid::set_global_grid(SGTraits::idx_t nx, SGTraits::idx_t ny, SGTraits::idx_t nz)
{
    global_.size[0] = nx;
    global_.size[1] = ny;
    global_.size[2] = nz;

    global_.start[0] = 0;
    global_.start[1] = 0;
    global_.start[2] = 0;

    global_.end[0] = nx;
    global_.end[1] = ny;
    global_.end[2] = nz;

    gridInitialized_ = true;

    local_ = sgix::grow_box(global_, nghost_[0], nghost_[1], nghost_[2]);
}

void StructGrid::set_partitions(int px, int py, int pz)
{
  const auto &pinfo = get_mpi();
  const int rank = pinfo.rank();
  const int nblocks = px * py * pz;

  if (nblocks != pinfo.size())
      throw std::runtime_error("Invalid partition specified for StructGrid");

  partitions_[0] = px;
  partitions_[1] = py;
  partitions_[2] = pz;

  // Determine local block information
  int p, q, r;

  p = (rank % px);
  q = ((rank - p) / px) % py;
  r = (rank - p - px * q) / (px * py);

  partindex_[0] = p;
  partindex_[1] = q;
  partindex_[2] = r;

  const auto *gsize = global_.size;

  for (int ii = 0; ii < StructBox::ndim; ++ii) {
    int rem = (gsize[ii] % partitions_[ii]);
    int quot = (gsize[ii] / partitions_[ii]);

    int size = quot + ((partindex_[ii] < rem) ? 1 : 0);
    int start = (partindex_[ii] * quot) + std::min(partindex_[ii], rem);
    int end = start + size;

    local_.nghost[ii] = nghost_[ii];
    local_.size[ii] = size + 2 * nghost_[ii];
    local_.start[ii] = start - nghost_[ii];
    local_.end[ii] = end + nghost_[ii];
  }

  partInitialized_ = true;
}

}  // nalu
}  // sierra
