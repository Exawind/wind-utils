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

#ifndef CHANNELFIELDS_H
#define CHANNELFIELDS_H

#include "PreProcessingTask.h"

namespace sierra {
namespace nalu {

/**
 * Initialize velocity fields for channel flow simulations
 *
 * This task is activated by using the `init_channel_fields` task in the
 * preprocessing input file. It requires a section `init_channel_fields` in the
 * `nalu_preprocess` section with the following parameters:
 *
 *  ```
 *  init_channel_fields:
 *    fluid_parts: [Unspecified-2-HEX]
 *
 *    velocity:
 *      Re_tau : 550
 *      viscosity : 0.0000157
 *  ```
 *
 * The user specified the friction Reynolds number, `Re_tau`, and the
 * kinematic `viscosity` (in m^2/s). The velocity field is initialized
 * to a Reichardt function, with an imposed sinusoidal perturbation
 * and random perturbation in the wall parallel directions.
 */
class ChannelFields: public PreProcessingTask
{
public:
    template<typename T>
    using Array2D = std::vector<std::vector<T>>;

    ChannelFields(CFDMesh&, const YAML::Node&);

    virtual ~ChannelFields() {}

    //! Declare velocity fields and register them for output
    void initialize();

    //! Initialize the velocity fields by linear interpolation
    void run();

private:
    ChannelFields() = delete;
    ChannelFields(const ChannelFields&) = delete;

    //! Parse the YAML file and initialize parameters
    void load(const YAML::Node&);

    //! Helper function to parse and initialize velocity inputs
    void load_velocity_info(const YAML::Node&);

    //! Helper function to parse and initialize tke inputs
    void load_tke_info(const YAML::Node&);

    //! Helper function to setup the channel flow parameters
    void setup_parameters();

    //! Reichardt function
    double reichardt(const double y);

    //! Perturbation function for u
    double u_perturbation(const double x, const double y, const double z);

    //! Perturbation function for w
    double w_perturbation(const double x, const double y, const double z);

    //! Initialize the velocity field through linear interpolation
    void init_velocity_field();

    //! Initialize the tke field through linear interpolation
    void init_tke_field();

    //! STK Metadata object
    stk::mesh::MetaData& meta_;

    //! STK Bulkdata object
    stk::mesh::BulkData& bulk_;

    //! Parts of the fluid mesh where velocity is initialized
    stk::mesh::PartVector fluid_parts_;

    //! Dimensionality of the mesh
    int ndim_;

    //! Flag indicating whether velocity is initialized
    bool doVelocity_;

    //! Flag indicating whether turbulent kinetic energy is initialized
    bool doTKE_;

    //! Bounding box length (x-direction)
    double length_;

    //! Bounding box length (z-direction)
    double width_;

    //! Bounding box length (y-direction)
    double height_;

    //! Reichardt function integration parameter
    double C_;

    //! Skin friction Reynolds number
    double Re_tau_;

    //! Skin friction velocity
    double utau_;

    //! Kinematic viscosity
    double viscosity_;

    //! Von Karman constant
    const double kappa_ = 0.4;

    //! Seed for RNG
    const int seed_ = 2864;

    //! Wavenumber of sinusoidal perturbation for u
    const int k_pert_u_ = 16;

    //! Wavenumber of sinusoidal perturbation for w
    const int k_pert_w_ = 16;

    //! Amplitude of sinusoidal perturbation for u
    const int a_pert_u_ = 10;

    //! Amplitude of sinusoidal perturbation for w
    const int a_pert_w_ = 10;

    //! Amplitude of random perturbation for u
    const double a_rand_u_ = 5;

    //! Amplitude of random perturbation for w
    const double a_rand_w_ = 5;
};

}
}

#endif /* CHANNELFIELDS_H */
