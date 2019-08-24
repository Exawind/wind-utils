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

#include "netcdf.h"

#include <iostream>

namespace sierra {
namespace nalu {

namespace {

/** Utility function to check NetCDF errors and generate appropriate runtime
 *  error messages.
 */
inline void check_nc_error(int ierr)
{
    if (ierr != NC_NOERR)
        throw std::runtime_error("NetCDF Error: " + std::string(nc_strerror(ierr)));
}

}

DEFINE_INHERITANCE_REGISTRY(TurbulenceFile);

TurbulenceFile*
TurbulenceFile::create(std::string fmt)
{
    auto it = TurbulenceFileReg_ConstructorTable_->find(fmt);

    if (it != TurbulenceFileReg_ConstructorTable_->end()) {
        return (it->second)();
    } else {
        std::cout << "ERROR: Invalid Turbulence file format = " << fmt << std::endl;
        std::cout << "Valid formats are: " << std::endl;
        for (const auto& t: *TurbulenceFileReg_ConstructorTable_) {
            std::cout << "\t" << t.first << std::endl;
        }
    }

    return nullptr;
}

void TurbulenceFile::write_netcdf(std::string outfile)
{
    int ierr;
    int ncid;

    int sDim, xDim, yDim, zDim;
    int uid, vid, wid;
    int urefid, hrefid, hoffsetid;
    int dxid, scaleid, profid;
    int profile_typ;

    const std::string title = this->title();
    const std::string velUnits = "m/s";
    const std::string lenUnits = "m";
    const std::string windFlagComment = "1 => Ux = Umean (z) + uf; 0 => No, (U = uf)";

    std::cout << "Begin output in NetCDF format: " << outfile << std::endl;
    ierr = nc_create(outfile.c_str(), NC_CLOBBER | NC_NETCDF4, &ncid);
    check_nc_error(ierr);

    // Wind profile type information
    {
        ierr = nc_def_enum(ncid, NC_BYTE, "wind_profile_t", &profile_typ);
        check_nc_error(ierr);

        signed char eval;

        eval = CONSTANT;
        ierr = nc_insert_enum(ncid, profile_typ, "Constant", &eval);
        eval = POWER_LAW;
        ierr = nc_insert_enum(ncid, profile_typ, "Power_Law", &eval);
        eval = LOG_LAW;
        ierr = nc_insert_enum(ncid, profile_typ, "Log_Law", &eval);
    }

    // Define dimensions for NetCDF file
    ierr = nc_def_dim(ncid, "ndim", static_cast<int>(ndim), &sDim);
    ierr = nc_def_dim(ncid, "nx", nx_[0], &xDim);
    ierr = nc_def_dim(ncid, "ny", nx_[1], &yDim);
    ierr = nc_def_dim(ncid, "nz", nx_[2], &zDim);

    // Type of wind profile
    ierr = nc_def_var(ncid, "wind_profile", profile_typ, 0, nullptr, &profid);

    // Scalar data
    ierr = nc_def_var(ncid, "uref", NC_DOUBLE, 0, NULL, &urefid);
    ierr = nc_def_var(ncid, "href", NC_DOUBLE, 0, NULL, &hrefid);
    ierr = nc_def_var(ncid, "hoffset", NC_DOUBLE, 0, NULL, &hoffsetid);

    int shearid, z0id;
    switch (profileType_) {
    case POWER_LAW:
        ierr = nc_def_var(ncid, "shear_exponent", NC_DOUBLE, 0, NULL, &shearid);
        break;

    case LOG_LAW:
        ierr = nc_def_var(ncid, "roughness_height", NC_DOUBLE, 0, NULL, &z0id);
        ierr = nc_put_att_text(ncid, z0id, "units", lenUnits.size() + 1, lenUnits.c_str());
        break;

    default:
        break;
    }

    int windFlag;
    ierr = nc_def_var(ncid, "mean_wind_added", NC_INT, 0, NULL, &windFlag);

    // Grid resolutions
    ierr = nc_def_var(ncid, "box_lengths", NC_DOUBLE, 1, &sDim, &dxid);

    ierr = nc_def_var(ncid, "scale_factors", NC_DOUBLE, 1, &sDim, &scaleid);

    // Define velocity arrays
    const std::vector<int> threeD { xDim, yDim, zDim };
    ierr = nc_def_var(ncid, "uvel", NC_DOUBLE, static_cast<int>(ndim),
                      threeD.data(), &uid);
    ierr = nc_def_var(ncid, "vvel", NC_DOUBLE, static_cast<int>(ndim),
                      threeD.data(), &vid);
    ierr = nc_def_var(ncid, "wvel", NC_DOUBLE, static_cast<int>(ndim),
                      threeD.data(), &wid);

    // Attributes
    ierr = nc_put_att_text(ncid, NC_GLOBAL, "title", title.size() + 1, title.c_str());

    ierr = nc_put_att_text(ncid, hrefid, "units", lenUnits.size() + 1, lenUnits.c_str());
    ierr = nc_put_att_text(ncid, hoffsetid, "units", lenUnits.size() + 1,
                           lenUnits.c_str());
    ierr = nc_put_att_text(ncid, dxid, "units", lenUnits.size() + 1, lenUnits.c_str());

    ierr = nc_put_att_text(ncid, urefid, "units", velUnits.size() + 1, velUnits.c_str());
    ierr = nc_put_att_text(ncid, uid, "units", velUnits.size() + 1, velUnits.c_str());
    ierr = nc_put_att_text(ncid, vid, "units", velUnits.size() + 1, velUnits.c_str());
    ierr = nc_put_att_text(ncid, wid, "units", velUnits.size() + 1, velUnits.c_str());

    ierr = nc_put_att_text(ncid, windFlag, "comment", windFlagComment.size() + 1,
                           windFlagComment.c_str());

    // Indicate that we are done defining variables and are ready to write data
    ierr = nc_enddef(ncid);
    check_nc_error(ierr);

    // Populate the data
    unsigned char proftype = profileType_;
    ierr = nc_put_var(ncid, profid, &proftype);

    switch (profileType_) {
    case POWER_LAW:
        ierr = nc_put_var(ncid, shearid, &shear_exp_);
        break;

    case LOG_LAW:
        ierr = nc_put_var(ncid, z0id, &z0_);
        break;

    default:
        break;
    }

    int meanWindAdded = add_mean_wind_ ? 1 : 0;

    ierr = nc_put_var(ncid, windFlag, &meanWindAdded);
    ierr = nc_put_var_double(ncid, dxid, boxlen_);
    ierr = nc_put_var_double(ncid, scaleid, scale_factors_);
    ierr = nc_put_var_double(ncid, urefid, &uref_);
    ierr = nc_put_var_double(ncid, hrefid, &href_);
    ierr = nc_put_var_double(ncid, hoffsetid, &hoffset_);
    ierr = nc_put_var_double(ncid, uid, uvel_.data());
    ierr = nc_put_var_double(ncid, vid, vvel_.data());
    ierr = nc_put_var_double(ncid, wid, wvel_.data());

    ierr = nc_close(ncid);
    check_nc_error(ierr);
    std::cout << "NetCDF file written successfully: " << outfile << std::endl;
}

}  // nalu
}  // sierra
