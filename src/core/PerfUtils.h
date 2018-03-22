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

#ifndef PERFUTILS_H
#define PERFUTILS_H

#include "Teuchos_TimeMonitor.hpp"
#include "stk_util/environment/perf_util.hpp"

namespace sierra {
namespace nalu {

/** Return a timer identified by name.
 *
 *  If an existing timer is found, then the timer is returned. Otherwise a new
 *  timer is created. The user will have to manually start/stop the timer. For
 *  most use cases, it might be preferable to use ``get_stopwatch`` function
 *  instead.
 */
inline Teuchos::RCP<Teuchos::Time> get_timer(const std::string& name)
{
    Teuchos::RCP<Teuchos::Time> timer = Teuchos::TimeMonitor::lookupCounter(name);
    if (timer.is_null())
        timer = Teuchos::TimeMonitor::getNewCounter(name);

    return timer;
}

/** Return a stopwatch identified by name.
 *
 *  The clock starts automatically upon invocation and will be stopped once the
 *  ``Teuchos::Timemonitor`` instance returned by this function goes out of
 *  scope.
 *
 */
inline Teuchos::TimeMonitor get_stopwatch(const std::string& name)
{
    return Teuchos::TimeMonitor(*get_timer(name));
}

std::string human_bytes(const size_t bytes);

void summarize_memory_usage(MPI_Comm comm, std::ostream& out = std::cout);

}  // nalu
}  // sierra


#endif /* PERFUTILS_H */
