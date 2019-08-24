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

#ifndef WINDSIMFILE_H
#define WINDSIMFILE_H

#include "tools/turbsim_netcdf/TurbulenceFile.h"

namespace sierra {
namespace nalu {

/** Wind Simulator (Mann turbulence model) interface
 *
 *  This class provides an interface to interact with turbulent fields generated
 *  from DTU's WindSim utility (Mann wind fields).
 */
class WindSimFile : public TurbulenceFile
{
public:
    WindSimFile() = default;

    virtual ~WindSimFile() = default;

    //! Load turbulence data from binary files
    virtual void load_turbulence_data(const YAML::Node&) override;

    //! Unique string identifier for this model
    virtual std::string title() override
    { return title_; }

private:
    //! Unique string identifier for the turbulence model
    std::string title_{"WindSim: Mann turbulence file"};
};

}  // nalu
}  // sierra


#endif /* WINDSIMFILE_H */
