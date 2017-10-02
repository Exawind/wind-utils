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

#ifndef CONSTANTSPACING_H
#define CONSTANTSPACING_H

#include "MeshSpacing.h"

namespace sierra {
namespace nalu {

/** Constant mesh spacing distribution
 *
 *  Specialization of MeshSpacing to allow for constant mesh spacing which is
 *  the default implementation if no user option is specified in the input file.
 *  This class requires no additional input arguments in the YAML file.
 */
class ConstantSpacing : public MeshSpacing
{
public:
    ConstantSpacing(
        int,
        const YAML::Node&);

    //! Initialize a constant spacing 1-D mesh
    virtual void init_spacings();

private:
    ConstantSpacing() = delete;
    ConstantSpacing(const ConstantSpacing&) = delete;
};

}  // nalu
}  // sierra



#endif /* CONSTANTSPACING_H */
