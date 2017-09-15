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

/** \file
 * Miscellaneous utilities for working with YAML C++ library
 */

#include "yaml-cpp/yaml.h"

namespace sierra {
namespace nalu {
namespace wind_utils {

template<typename T>
bool get_optional(const YAML::Node& node, const std::string& key, T& result)
{
    bool found = false;
    if (node[key]) {
        result = node[key].as<T>();
        found = true;
    }

    return found;
}

template<typename T>
bool get_optional(const YAML::Node& node, const std::string& key, T& result, const T& default_value)
{
    auto found = get_optional(node, key, result);
    if (!found) {
        result = default_value;
    }
    return found;
}

}  // wind_utils
}  // nalu
}  // sierra
