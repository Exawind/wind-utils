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

#include "ConstantSpacing.h"

namespace sierra {
namespace nalu {

REGISTER_DERIVED_CLASS(MeshSpacing, ConstantSpacing, "constant_spacing");

ConstantSpacing::ConstantSpacing(
    int npts,
    const YAML::Node&
) : MeshSpacing(npts)
{
    // Ignore yaml node... no additional inputs
}

void ConstantSpacing::init_spacings()
{
    double dx = 1.0 / static_cast<double>(numPts_-1);

    for (int i=0; i<numPts_; i++)
        ratios_[i] = i * dx;
}

}  // nalu
}  // sierra
