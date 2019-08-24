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

#ifndef TURBULENCEFILE_H
#define TURBULENCEFILE_H

#include "core/ClassRegistry.h"
#include "yaml-cpp/yaml.h"

#include <string>
#include <vector>

namespace sierra {
namespace nalu {

/** Synthetic turbulence interface
 *
 *  This class represents a generic interface to synthetic turbulence data
 *  (e.g., Mann wind fields, TurbSim, HIT, etc.) that can be used with ExaWind
 *  framework.
 */
class TurbulenceFile
{
public:
    //! Spatial dimension
    static constexpr unsigned int ndim = 3;

    //! Type of wind profile that is superimposed on the turbulent field
    enum WindProfileType {
        CONSTANT = 0, ///< Constant wind speed
        POWER_LAW,    ///< Power-law `U = U_\inf (z/H)^\alpha`
        LOG_LAW       ///< Log-law
    };

    TurbulenceFile() = default;

    virtual ~TurbulenceFile() = default;

    //! Load turbulence from files
    virtual void load_turbulence_data(const YAML::Node&) = 0;

    //! A string representation of the type of turbulence field
    virtual std::string title() = 0;

    //! Output the turbulence file in NetCDF format
    void write_netcdf(std::string);

    DECLARE_INHERITANCE_REGISTRY
    (
        TurbulenceFile,
        (),
        ()
    );

    /** Create and return a concrete instance of TurbulenceFile reader
     *
     *  @param fileType The format of the turbulence file
     *  @return Instance of a turbulence file reader
     */
    static TurbulenceFile* create(std::string);

protected:
    //! u component [nx, ny, nz]
    std::vector<double> uvel_;

    //! v component [nx, ny, nz]
    std::vector<double> vvel_;

    //! w component [nx, ny, nz]
    std::vector<double> wvel_;

    //! Box length in each direction
    double boxlen_[ndim];

    //! Grid size in each direction
    double dx_[ndim];

    //! Turbulence scaling factors
    double scale_factors_[ndim] = {1.0, 1.0, 1.0};

    //! Number of grid points in each direction
    int nx_[ndim];

    //! Reference velocity (m/s)
    double uref_;

    //! Reference height (m)
    double href_;

    //! Offset height
    double hoffset_{0.0};

    //! Shear exponent for power-law profiles
    double shear_exp_{0.0};

    //! Roughness height for log-law profiles (m)
    double z0_{0.1};

    //! Total number of grid points
    unsigned int npts_{1};

    //! Wind profile type
    WindProfileType profileType_{CONSTANT};

    //! Flag indicating whether a mean wind profile has been added to the field
    bool add_mean_wind_{true};
};

}  // nalu
}  // sierra


#endif /* TURBULENCEFILE_H */
