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

#ifndef SLICE_H
#define SLICE_H

#include "core/CFDMesh.h"
#include "core/YamlUtils.h"

#include "Kokkos_Macros.hpp"
#include "Kokkos_Core.hpp"

namespace sierra {
namespace nalu {

class Slice
{
public:
    static constexpr int NDim = 3;
    using AxesType =
        Kokkos::View<double[NDim][NDim], Kokkos::DefaultHostExecutionSpace>;
    using BBoxType =
        Kokkos::View<double[4][NDim], Kokkos::DefaultHostExecutionSpace>;

    Slice(CFDMesh&);

    virtual ~Slice() = default;

    virtual void load(const YAML::Node&);

    virtual void initialize();

    virtual void run();

protected:
    CFDMesh& mesh_;

    AxesType axes_{"axes"};

    BBoxType vertices_{"vertices"};

    std::vector<double> grid_dx_;

    std::vector<size_t> grid_dims_;

    std::vector<double> plane_offsets_;

    stk::mesh::PartVector partVec_;

    std::string partNamePrefix_{"plane_"};

    unsigned num_planes_{1};
};

}  // nalu
}  // sierra


#endif /* SLICE_H */
