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

#ifndef ABLFIELDS_H
#define ABLFIELDS_H

#include "PreProcessingTask.h"

namespace sierra {
namespace nalu {

/**
 * Initialize velocity and temperature fields for ABL simulations
 *
 * This task is activated by using the `init_abl_fields` task in the
 * preprocessing input file. It requires a section `init_abl_fields` in the
 * `nalu_preprocess` section with the following parameters:
 *
 *  ```
 *  init_abl_fields:
 *    fluid_parts: [Unspecified-2-HEX]
 *
 *    temperature:
 *      heights: [    0, 650.0, 750.0, 10750.0]
 *      values:  [280.0, 280.0, 288.0,   318.0]
 *
 *    velocity:
 *      heights: [0.0, 10.0, 30.0, 70.0, 100.0, 650.0, 10000.0]
 *      values:
 *        - [ 0.0, 0.0, 0.0]
 *        - [4.81947, -4.81947, 0.0]
 *        - [5.63845, -5.63845, 0.0]
 *        - [6.36396, -6.36396, 0.0]
 *        - [6.69663, -6.69663, 0.0]
 *        - [8.74957, -8.74957, 0.0]
 *        - [8.74957, -8.74957, 0.0]
 *  ```
 *
 * The sections `temperature` and `velocity` are optional, allowing the user to
 * initialize only the temperature or the velocity as desired. The heights are
 * in meters, the temperature is the potential temperature in Kelvin, and the
 * velocity is the actual vector in m/s. Currently, the code does not include
 * the ability to automatically convert (mangitude, direction) to velocity
 * vectors.
 */
class ABLFields: public PreProcessingTask
{
public:
    template<typename T>
    using Array2D = std::vector<std::vector<T>>;

    /**
     * \param mesh A sierra::nalu::CFDMesh instance
     * \param node The YAML::Node containing inputs for this task
     */
    ABLFields(CFDMesh&, const YAML::Node&);

    virtual ~ABLFields() {}

    //! Declare velocity and temperature fields and register them for output
    void initialize();

    //! Initialize the velocity and/or temperature fields by linear interpolation
    void run();

private:
    ABLFields() = delete;
    ABLFields(const ABLFields&) = delete;

    //! Parse the YAML file and initialize parameters
    void load(const YAML::Node&);

    //! Helper function to parse and initialize velocity inputs
    void load_velocity_info(const YAML::Node&);

    //! Helper function to parse and initialize temperature inputs
    void load_temperature_info(const YAML::Node&);

    //! Initialize the velocity field through linear interpolation
    void init_velocity_field();

    //! Intialize the temperature field through linear interpolation
    void init_temperature_field();

    //! Add perturbations to velocity field
    void perturb_velocity_field();

    //! Add perturbations to temperature field
    void perturb_temperature_field();

    //! STK Metadata object
    stk::mesh::MetaData& meta_;

    //! STK Bulkdata object
    stk::mesh::BulkData& bulk_;

    //! Parts of the fluid mesh where velocity/temperature is initialized
    stk::mesh::PartVector fluid_parts_;

    //! List of heights where velocity is defined
    std::vector<double> vHeights_;

    //! List of velocity (3-d components) at the user-defined heights
    Array2D<double> velocity_;

    //! List of heights where temperature is defined
    std::vector<double> THeights_;

    //! List of temperatures (K) at user-defined heights (THeights_)
    std::vector<double> TValues_;

    //! List of periodic parts
    std::vector<std::string> periodicParts_;

    //! Velocity perturbation amplitude for Ux
    double deltaU_{1.0};

    //! Velocity perturbation amplitude for Uy
    double deltaV_{1.0};

    //! Number of periods for Ux
    double Uperiods_{4.0};

    //! Number of periods for Uy
    double Vperiods_{4.0};

    //! Reference height for velocity perturbations
    double zRefHeight_{50.0};

    //! Amplitude of temperature perturbations
    double thetaAmplitude_;

    //! Mean for the Gaussian random number generator
    double thetaGaussMean_{0.0};

    //! Variance of the Gaussian random number generator
    double thetaGaussVar_{1.0};

    //! Cutoff height for temperature fluctuations
    double thetaCutoffHt_;

    //! Dimensionality of the mesh
    int ndim_;

    //! Flag indicating whether velocity is initialized
    bool doVelocity_;

    //! Flag indicating whether temperature is initialized
    bool doTemperature_;

    //! Flag indicating whether velocity perturbations are added during initialization
    bool perturbU_{false};

    //! Flag indicating whether temperature perturbations are added
    bool perturbT_{false};
};

}
}

#endif /* ABLFIELDS_H */
