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

#include "GeometricStretching.h"
#include "core/YamlUtils.h"

#include <iostream>
#include <cmath>

namespace sierra {
namespace nalu {

REGISTER_DERIVED_CLASS(MeshSpacing, GeometricStretching, "geometric_stretching");

GeometricStretching::GeometricStretching(
    int npts,
    const YAML::Node& node
) : MeshSpacing(npts)
{
    load(node);
}

void GeometricStretching::load(const YAML::Node& node)
{
    wind_utils::get_optional(node, "bidirectional", bidirectional_);
    wind_utils::get_optional(node, "verbose", verbose_);
    fac_ = node["stretching_factor"].as<double>();

    int ny = bidirectional_? (numPts_ - 1)/2 + (numPts_ - 1) % 2: (numPts_ - 1);
    double tlen = bidirectional_? 0.5 : 1.0;
    fch_ = tlen * (fac_ - 1.0) / (std::pow(fac_, ny) - 1.0);
}

void GeometricStretching::init_spacings()
{
    if (!bidirectional_)
        unidirectional_stretching();
    else
        bidirectional_stretching();

    if (verbose_)
        for (int i=0; i < numPts_; i++)
            std::cerr << i << "\t" << ratios_[i] << std::endl;
}

void GeometricStretching::bidirectional_stretching()
{
    ratios_[0] = 0.0;
    ratios_[numPts_ - 1] = 1.0;
    int ny = (numPts_ - 1) / 2;

    double rfac = 1.0;
    for (int i=1; i <= ny; i++) {
        ratios_[i] = ratios_[i-1] + fch_ * rfac;
        ratios_[numPts_ - 1 - i] = ratios_[numPts_ - i] - fch_ * rfac;
        rfac *= fac_;
    }

}

void GeometricStretching::unidirectional_stretching()
{
    ratios_[0] = 0;
    double rfac = 1.0;
    for (int i=1; i < numPts_; i++) {
        ratios_[i] = ratios_[i-1] + fch_ * rfac;
        rfac *= fac_;
    }
}

}  // nalu
}  // sierra
