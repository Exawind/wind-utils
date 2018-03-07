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

/**\file
 * Utilities to use Kokkos loops within NaluWindUtils
 */
#ifndef KOKKOSWRAPPERS_H
#define KOKKOSWRAPPERS_H

#include "Kokkos_Core.hpp"

namespace sierra {
namespace nalu {

using DynamicScheduleType = Kokkos::Schedule<Kokkos::Dynamic>;
using TeamPolicyType = Kokkos::TeamPolicy<>;
using RangePolicyType = Kokkos::RangePolicy<DynamicScheduleType>;
using TeamMemberType = Kokkos::TeamPolicy<DynamicScheduleType>::member_type;

}  // nalu
}  // sierra


#endif /* KOKKOSWRAPPERS_H */
