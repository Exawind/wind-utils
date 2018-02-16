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

#ifndef NESTEDREFINEMENT_H
#define NESTEDREFINEMENT_H

#include "PreProcessingTask.h"

#include <vector>
#include <array>

namespace sierra {
namespace nalu {

/** Tag regions in mesh for refinement with Percept mesh_adapt utility.
 *
 *  This utility creates a field turbine_refinement_field that is populated with
 *  an indicator value between [0, 1] that can be used with the Percept
 *  mesh_adapt utility to locally refine regions of interest.
 *
 *  A typical use of this utility is to refine an ABL mesh around turbines,
 *  especially for use with actuator line wind farm simulations.
 */
class NestedRefinement: public PreProcessingTask
{
public:
    typedef std::array<double, 3> Vec3D;
    typedef std::array<std::array<double, 3>, 3> TrMat;

    NestedRefinement(CFDMesh&, const YAML::Node&);

    virtual ~NestedRefinement() {}

    //! Initialize the refinement field and register to parts
    void initialize();

    //! Perform search and tag elements with appropriate values for subsequent
    //! refinement
    void run();

private:
    NestedRefinement() = delete;
    NestedRefinement(const NestedRefinement&) = delete;

    //! Parse the YAML file and initialize the necessary parameters
    void load(const YAML::Node&);

    //! Process input data and populate necessary data structures for subsequent use
    void process_inputs();

    //! Estimate the refinement fraction [0,1] for a given element, indicated by
    //! the element mid point.
    double compute_refine_fraction(Vec3D&);

    //! Write out the input files that can be used with Percept
    void write_percept_inputs();

    //! Partnames for the ABL mesh
    std::vector<std::string> fluidPartNames_;

    //! Parts of the ABL mesh where refinement is performed
    stk::mesh::PartVector fluidParts_;

    //! List of turbine diameters for the turbines in the wind farm [numTurbines]
    std::vector<double> turbineDia_;

    //! List of turbine tower heights for the turbines in wind farm [numTurbines]
    std::vector<double> turbineHt_;

    //! List of turbine pad locations [numTurbines, 3]
    std::vector<Vec3D> turbineLocs_;

    //! List of refinement levels [numLevels, 3]
    std::vector<std::vector<double>> refineLevels_;

    //! Transformation matrices for each turbine [numTurbines]
    std::vector<TrMat> boxAxes_;

    //! The minimum corners for each refinement box [numTurbines * numLevels]
    std::vector<Vec3D> corners_;

    //! The dimensions of each box [numTurbines * numLevels]
    std::vector<Vec3D> boxLengths_;

    //! Field name used in the Exodus mesh for the error indicator field
    std::string refineFieldName_{"turbine_refinement_field"};

    //! Prefix for the input file name
    std::string perceptFilePrefix_{"adapt"};

    //! Search tolerance used when searching for box inclusion
    double searchTol_{10.0};

    //! Compass direction of the wind (in degrees)
    double windAngle_{270.0};

    //! The number of turbines in the wind farm
    size_t numTurbines_;

    //! The number of refinement levels
    size_t numLevels_;

    //! Write input files for use with subsequent percept run
    bool writePercept_{true};
};

}  // nalu
}  // sierra


#endif /* NESTEDREFINEMENT_H */
