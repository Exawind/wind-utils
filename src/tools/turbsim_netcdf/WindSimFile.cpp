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

#include "tools/turbsim_netcdf/WindSimFile.h"
#include "core/YamlUtils.h"

#include <cassert>
#include <fstream>
#include <algorithm>
#include <iostream>

namespace sierra {
namespace nalu {

REGISTER_DERIVED_CLASS(TurbulenceFile, WindSimFile, "windsim");

namespace {

/** Load a WindSim binary file
 *
 *  @param fname Binary filename
 *  @param buffer Float (4-byte real) array to read in the data
 *  @param nbytes Total number of bytes to be read
 */
inline void load_bin_file(std::string fname, std::vector<float>& buffer, unsigned int nbytes)
{
    std::ifstream binfile(fname, std::ios::in | std::ios::binary);

    if (!binfile.is_open())
        throw std::runtime_error("WindSimFile: Error opening file = " + fname);
    else
        std::cout << "\tLoading file: " << fname << std::endl;

    binfile.read(reinterpret_cast<char*>(buffer.data()), nbytes);

    binfile.close();
}

}

void WindSimFile::load_turbulence_data(const YAML::Node &node)
{
    std::cout << "Begin loading WindSim turbulence data" << std::endl;
    auto box_dims = node["box_dims"].as<std::vector<int>>();
    auto box_len = node["box_len"].as<std::vector<double>>();
    auto filenames = node["bin_filenames"].as<std::vector<std::string>>();
    assert(box_dims.size() == ndim);
    assert(box_len.size() == ndim);
    assert(filenames.size() == ndim);

    {
      const auto &wnode = node["wind_profile"];
      uref_ = wnode["ref_speed"].as<double>();
      href_ = wnode["ref_height"].as<double>();

      wind_utils::get_optional(wnode, "height_offset", hoffset_);

      if (wnode["profile_type"]) {
          const std::string proftype = wnode["profile_type"].as<std::string>();

          if (proftype == "power_law")
              profileType_ = POWER_LAW;
          else if (proftype == "log_law")
              profileType_ = LOG_LAW;
          else if (proftype == "constant")
              profileType_ = CONSTANT;
          else
              throw std::runtime_error("Invalid profile type: " + proftype);

          switch (profileType_) {
          case POWER_LAW:
              shear_exp_ = wnode["shear_exponent"].as<double>();
              break;

          case LOG_LAW:
              z0_ = wnode["roughness_height"].as<double>();

          default:
              break;
          }
      }
    }

    npts_ = 1;
    for (unsigned int i=0; i < ndim; ++i) {
        nx_[i] = box_dims[i];
        boxlen_[i] = box_len[i];

        // From wind sim documentation, the distance between two points is
        // (L_i/N_i) where i=1..3
        dx_[i] = box_len[i] / static_cast<double>(nx_[i]);

        npts_ *= nx_[i];
    }
    uvel_.resize(npts_);
    vvel_.resize(npts_);
    wvel_.resize(npts_);

    const auto nbytes = npts_ * sizeof(float);
    std::vector<float> buffer(npts_);

    // U component
    load_bin_file(filenames[0], buffer, nbytes);
    std::copy(buffer.begin(), buffer.end(), uvel_.begin());

    // V component
    load_bin_file(filenames[1], buffer, nbytes);
    std::copy(buffer.begin(), buffer.end(), vvel_.begin());

    // W component
    load_bin_file(filenames[2], buffer, nbytes);
    std::copy(buffer.begin(), buffer.end(), wvel_.begin());

    // Deal with scaling factors
    std::vector<double> scale_factors(ndim, 1.0);
    bool apply_scaling = false;

    std::string scale_type = "none";
    wind_utils::get_optional(node, "scale_type", scale_type);

    // Get scaling factors
    if (scale_type != "none") {
      scale_factors = node["scaling_factors"].as<std::vector<double>>();
      std::copy(scale_factors.begin(), scale_factors.end(), std::begin(scale_factors_));

      // Assume user wants to apply scaling
      apply_scaling = true;
    }

    wind_utils::get_optional(node, "apply_scaling", apply_scaling);
    if (apply_scaling) {
      std::cout << "\tApplying scale factors" << std::endl;
      for (unsigned int i = 0; i < npts_; i++) {
        uvel_[i] *= scale_factors[0];
        vvel_[i] *= scale_factors[1];
        wvel_[i] *= scale_factors[2];
      }
    }

    wind_utils::get_optional(node, "add_mean_wind", add_mean_wind_);

    if (add_mean_wind_) {
        std::cout << "Adding mean wind profile" << std::endl;
        auto zbase = href_ - 0.5 * boxlen_[2];
        auto dh = boxlen_[2] / static_cast<double>(nx_[2] - 1);

        double uvel = 0.0;
        for (int k = 0; k < nx_[2]; ++k) {
            const double ht = zbase + k * dh;

            switch (profileType_) {
            case POWER_LAW:
                uvel = uref_ * std::pow((ht - hoffset_) / href_, shear_exp_);
                break;

            case LOG_LAW:
                uvel = uref_ * (std::log((ht - hoffset_) / z0_) /
                                std::log(href_ / z0_));
                break;

            default:
                uvel = uref_;
                break;
            }

            for (int i=0; i < nx_[0]; ++i) {
                const int offset = i * nx_[1] * nx_[2];
                for (int j = 0; j < nx_[1]; ++ j) {
                    const int idx = offset + j * nx_[2];

                    uvel_[idx + k] += uvel;
                }
            }
        }
    }
}

}  // nalu
}  // sierra
