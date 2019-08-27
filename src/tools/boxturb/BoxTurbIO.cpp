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

#include <cassert>
#include <fstream>

#include "tools/boxturb/BoxTurbIO.h"
#include "core/ParallelInfo.h"
#include "struct_grid/StructGridIx.h"

namespace sierra {
namespace nalu {

namespace {

/** Load a WindSim binary file
 *
 *  @param fname Binary filename
 *  @param buffer Float (4-byte real) array to read in the data
 *  @param nbytes Total number of bytes to be read
 */
inline void load_bin_file(
    std::string fname,
    std::vector<float> &buffer,
    unsigned int nbytes)
{
    std::ifstream binfile(fname, std::ios::in | std::ios::binary);

    if (!binfile.is_open())
        throw std::runtime_error("WindSimFile: Error opening file = " + fname);
    else
        get_mpi().info() << "\tLoading file: " << fname << std::endl;

    binfile.read(reinterpret_cast<char *>(buffer.data()), nbytes);

    binfile.close();
}
}

void
BoxTurbIO::load(const std::string& ftype, const YAML::Node& node)
{
    if (ftype == "windsim")
        load_windsim_file(node);
    else
        throw std::runtime_error("Invalid turbulence file format: " + ftype);
}

void
BoxTurbIO::load_windsim_file(const YAML::Node& node)
{
    using idx_t = SGTraits::idx_t;
    const auto& pinfo = get_mpi();

    pinfo.info()
        << "Begin loading WindSim turbulence data" << std::endl;

    auto filenames = node["bin_filenames"].as<std::vector<std::string>>();
    assert(filenames.size() == SGTraits::ndim);

    const auto ncells_global = sgix::num_cells(boxturb_.grid_.global());
    const auto nbytes = ncells_global * sizeof(float);
    std::vector<float> buffer(ncells_global);

    auto index = sgix::PeriodicIndexer<BoxTurb::Layout>(boxturb_.grid_.global());

    // U component
    {
        load_bin_file(filenames[0], buffer, nbytes);
        auto& uvel = *boxturb_.uvel_;
        SGTraits::size_t ii = 0;

        sgix::ijk_loop(
            boxturb_.grid_.local(),
            [&](idx_t i, idx_t j, idx_t k) {
                uvel[ii++] = buffer[index(i, j, k)];
            });
    }

    // V component
    {
        load_bin_file(filenames[1], buffer, nbytes);
        auto& vvel = *boxturb_.vvel_;
        SGTraits::size_t ii = 0;

        sgix::ijk_loop(
            boxturb_.grid_.local(),
            [&](idx_t i, idx_t j, idx_t k) {
                vvel[ii++] = buffer[index(i, j, k)];
            });
    }

    // W component
    {
        load_bin_file(filenames[2], buffer, nbytes);
        auto& wvel = *boxturb_.wvel_;

        SGTraits::size_t ii = 0;
        sgix::ijk_loop(
            boxturb_.grid_.local(),
            [&](idx_t i, idx_t j, idx_t k) {
                wvel[ii++] = buffer[index(i, j, k)];
            });
    }

    boxturb_.source_ = "WindSim Mann turbulence field";
}

}  // nalu
}  // sierra
