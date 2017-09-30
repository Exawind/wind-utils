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

#ifndef BDYIOPLANES_H
#define BDYIOPLANES_H

#include "PreProcessingTask.h"

namespace sierra {
namespace nalu {

/** Extract boundary planes for I/O mesh.
 *
 *  Given an ABL precursor mesh, this utility extracts the specified boundaries
 *  and creates a new IO Transfer mesh for use with ABL precursor simulations.
 *
 */
class BdyIOPlanes: public PreProcessingTask
{
public:
    /**
     * \param mesh A sierra::nalu::CFDMesh instance
     * \param node The YAML::Node containing inputs for this task
     */
    BdyIOPlanes(CFDMesh&, const YAML::Node&);

    virtual ~BdyIOPlanes();

    /** Register boundary parts and attach coordinates to the parts
     *
     *  The parts are created as SHELL elements to as needed by the Nalu
     *  Transfer class.
     */
    virtual void initialize();

    /** Copy user specified boundaries and save the IO Transfer mesh.
     */
    virtual void run();

private:
    BdyIOPlanes() = delete;
    BdyIOPlanes(const BdyIOPlanes&) = delete;

    //! Parse user inputs from the YAML file
    void load(const YAML::Node&);

    //! Copy the boundary from Fluid mesh to the IO Xfer mesh
    void create_boundary(const std::string);

    //! Original mesh DB information
    CFDMesh& mesh_;

    //! IO Mesh db STK meta and bulk data
    CFDMesh iomesh_;

    //! User specified list of boundaries to be extracted
    std::vector<std::string> bdyNames_;

    //! Name of the I/O db where the boundaries are written out
    std::string output_db_{""};
};

}  // nalu
}  // sierra

#endif /* BDYIOPLANES_H */
