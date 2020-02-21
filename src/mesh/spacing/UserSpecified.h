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

#ifndef USERSPECIFIED_H
#define USERSPECIFIED_H

#include "MeshSpacing.h"

namespace sierra {
namespace nalu {

/** Create a mesh spacing distribution from user inputs
 *
 *  Requires the user to specify an array of node positions over a unit span
 */
class UserSpecified : public MeshSpacing
{
public:
    UserSpecified() = delete;
    UserSpecified(const UserSpecified&) = delete;

    UserSpecified(int, const YAML::Node&);

    // Initialize the spacing (required method)
    virtual void init_spacings();

private:
    //! Process the user-provided YAML inputs
    void load(const YAML::Node&);
};

}  // nalu
}  // sierra


#endif /* USERSPECIFIED_H */
