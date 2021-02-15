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

#include "tools/turbsim_netcdf/TurbulenceFile.h"
#include "core/YamlUtils.h"

#include "stk_util/parallel/Parallel.hpp"
#include "stk_util/environment/OptionsSpecification.hpp"
#include "stk_util/environment/ParseCommandLineArgs.hpp"
#include "stk_util/environment/ParsedOptions.hpp"

#include "netcdf.h"
#include "netcdf_par.h"

#include <iostream>
#include <fstream>
#include <cctype>
#include <algorithm>
#include <memory>

int backup(int argc, char** argv)
{
    stk::ParallelMachine comm = stk::parallel_machine_init(&argc, &argv);

    int nproc = stk::parallel_machine_size(comm);
    if (nproc > 1)
        throw std::runtime_error("Cannot run turbulence file processor in parallel");

    std::string filename;
    stk::OptionsSpecification desc(
        "Nalu turbsim convertor utility. Valid options are");
    desc.add_options()("help,h", "Show this help message")
        ("input-file,i", "Input file for turbsim",
         stk::TargetPointer<std::string>(&filename),
         stk::DefaultValue<std::string>("turbsim_netcdf.yaml"));

    stk::ParsedOptions vmap;
    stk::parse_command_line_args(
        argc, const_cast<const char**>(argv), desc, vmap);

    if (vmap.count("help")) {
        std::cout << desc << std::endl;
        return 0;
    }

    filename = vmap["input-file"].as<std::string>();
    std::ifstream fin(filename.c_str());
    if (!fin.good()) {
        std::cout << "Cannot find input file: " << filename << std::endl;
        return 1;
    }

    std::cout << "\nNalu Turbulent File Processing Utility\n"
              << "Input file: " << filename << std::endl;

    YAML::Node inpfile(YAML::LoadFile(filename));
    const auto& node = inpfile["turbsim_netcdf"];
    const std::string output = node["output"].as<std::string>();

    std::string format = "windsim";
    sierra::nalu::wind_utils::get_optional(node, "data_format", format);

    std::transform(format.begin(), format.end(), format.begin(),
                   [](unsigned char c) { return std::tolower(c); });

    if (format == "hawc2") format = "windsim";
    if ((format == "fast") || (format == "openfast") || (format == "fast8"))
        format = "turbsim";

    using TurbulenceFile = sierra::nalu::TurbulenceFile;
    std::unique_ptr<TurbulenceFile> turbFile(TurbulenceFile::create(format));

    if (!turbFile)
        return 1;

    std::cout << "Turbulence file format: " << turbFile->title() << std::endl;
    turbFile->load_turbulence_data(node);
    turbFile->write_netcdf(output);

    stk::parallel_machine_finalize();
    return 0;
}

inline void
check_nc_error(int ierr)
{
    if (ierr != NC_NOERR)
        throw std::runtime_error(
            "NetCDF Error: " + std::string(nc_strerror(ierr)));
}

#define BTNC_CALL(funcarg)                      \
    do {                                        \
        int ierr = funcarg;                     \
        check_nc_error(ierr);                   \
    } while (0)


int main(int argc, char** argv)
{
    stk::ParallelMachine comm = stk::parallel_machine_init(&argc, &argv);

    // int iproc = stk::parallel_machine_size(comm);

    const std::string outfile = "test_ncwrite.nc";
    int ncid, ierr;
    ierr = nc_create_par(
        outfile.c_str(), NC_CLOBBER | NC_NETCDF4 | NC_MPIIO,
        comm, MPI_INFO_NULL, &ncid);
    check_nc_error(ierr);

    int recDim, timeID;
    ierr = nc_def_dim(ncid, "num_timesteps", NC_UNLIMITED, &recDim);
    ierr = nc_def_var(ncid, "time", NC_DOUBLE, 1, &recDim, &timeID);
    ierr = nc_var_par_access(ncid, timeID, NC_COLLECTIVE);
    ierr = nc_enddef(ncid);
    check_nc_error(ierr);

    size_t count0 = 1;
    for (size_t i=0; i < 10; ++i) {
        double time = 0.1 * i;
        ierr = nc_put_vara_double(ncid, timeID, &i, &count0, &time);
        check_nc_error(ierr);
    }

    ierr = nc_close(ncid);
    check_nc_error(ierr);
    stk::parallel_machine_finalize();
    return 0;
}
