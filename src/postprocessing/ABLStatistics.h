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

#ifndef ABLSTATISTICS_H
#define ABLSTATISTICS_H

#include "PostProcessingTask.h"

#include <unordered_map>

namespace sierra {
namespace nalu {

/** Process Atmospheric Boundary Layer statistics from precursor simulations
 *
 *  Extracts ABL profile information from a precursor run and calculates spatial
 *  averages of various quantities. The output is provided in tabular format in
 *  ASCII text files.
 *
 *  Current limitations:
 *    - Assumes constant spacing in z-direction
 */
class ABLStatistics: public PostProcessingTask
{
public:
    ABLStatistics(CFDMesh&, const YAML::Node&);

    virtual ~ABLStatistics() {}

    //! Declare necessary fields and initialize data structures for gathering statistics
    void initialize() override;

    //! Perform actions once the mesh is properly loaded
    void run() override;

private:
    ABLStatistics() = delete;
    ABLStatistics(const ABLStatistics&) = delete;

    //! Parse the YAML file and initialize the parameters
    void load(const YAML::Node&);

    /** Load solution fields for the final timestep available in the database
     *
     */
    void populate_solution();

    /** Compute and aggregate global spatial averages for quantities of interest
     */
    void average_planes();

    //! Output data
    void output_averages();

    //! STK Metadata object
    stk::mesh::MetaData& meta_;

    //! STK Bulkdata object
    stk::mesh::BulkData& bulk_;

    //! Parts of the fluid mesh where velocity/temperature is initialized
    stk::mesh::PartVector fluid_parts_;

    //! Mapping of required fields to actual field names
    std::unordered_map<std::string, std::string> field_map_;

    //! Spatially averaged mean velocity as a function of height [nHeights, 3]
    std::vector<double> velMean_;

    //! Spatially averaged mean temperature as a function of height [nHeights]
    std::vector<double> tempMean_;

    //! Spatially averaged mean SFS field as a function of height [nHeights, 6]
    std::vector<double> sfsMean_;

    //! Array of heights where averaged data is available
    std::vector<double> heights_;

    //! Running counter of nodes at a particular height for averaging purposes
    std::vector<int> node_counters_;

    //! Starting value of vertical height (of the domain). Default = 0.0
    double zmin_;

    //! Height of the top boundary of the ABL domain (must be specified by user)
    double zmax_;

    //! Constant spacing in the vertical direction for the ABL mesh
    double dz_;

    //! Total number of heights where data is sampled and averaged
    int nheights_;

    //! Dimensionality of the mesh (should be 3)
    const int ndim_{3};
};

}  // nalu
}  // sierra


#endif /* ABLSTATISTICS_H */
