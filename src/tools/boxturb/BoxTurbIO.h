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

#ifndef BOXTURBIO_H
#define BOXTURBIO_H

#include "tools/boxturb/BoxTurb.h"
#include "core/YamlUtils.h"

namespace sierra {
namespace nalu {

class BoxTurbIO
{
public:
    BoxTurbIO(BoxTurb& bt) : boxturb_(bt) {}

    virtual ~BoxTurbIO() = default;

    void load(const std::string&, const YAML::Node&);

    void load_windsim_file(const YAML::Node&);

protected:
    BoxTurb& boxturb_;
};

}  // nalu
}  // sierra


#endif /* BOXTURBIO_H */
