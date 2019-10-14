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

#ifndef CFDMESH_H
#define CFDMESH_H

#include "stk_util/parallel/Parallel.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/Entity.hpp"
#include "stk_mesh/base/Field.hpp"
#include "stk_mesh/base/CoordinateSystems.hpp"
#include "stk_search/Point.hpp"
#include "stk_search/Box.hpp"
#include "stk_io/StkMeshIoBroker.hpp"

#include <string>
#include <memory>
#include <unordered_set>

namespace sierra {
namespace nalu {

typedef stk::mesh::Field<double, stk::mesh::Cartesian> VectorFieldType;
typedef stk::mesh::Field<double> ScalarFieldType;
typedef stk::search::Point<double> PointType;
typedef stk::search::Box<double> BoxType;

/** STK Mesh interface
 *
 *  This class provides a thin wrapper around the STK mesh objects (MetaData,
 *  BulkData, and StkMeshIoBroker) for use with various preprocessing utilities.
 */
class CFDMesh
{
public:
    /**Create a CFD mesh instance from an existing mesh database
     *
     * \param comm MPI Communicator object
     * \param filename Exodus database filename
     */
    explicit CFDMesh(stk::ParallelMachine&,
            const std::string filename);

    /**Create a CFD mesh instance from scratch
     *
     * \param comm MPI Communicator object
     * \param ndim Dimensionality of mesh
     */
    explicit CFDMesh(stk::ParallelMachine&,
                     const int ndim);

    ~CFDMesh() {}

    /** Initialize the mesh database
     *
     *  If an input DB is provided, the mesh is read from the file. The MetaData
     *  is committed and the BulkData is ready for use/manipulation.
     */
    void init(stk::io::DatabasePurpose db_purpose = stk::io::READ_MESH);

    //! Reference to the MPI communicator object
    inline stk::ParallelMachine& comm() { return comm_; }

    //! Reference to the stk::mesh::MetaData instance
    inline stk::mesh::MetaData& meta() { return meta_; }

    //! Reference to the stk::mesh::BulkData instance
    inline stk::mesh::BulkData& bulk() { return bulk_; }

    //! Reference to the STK mesh I/O instance
    inline stk::io::StkMeshIoBroker& stkio() { return stkio_; }

    /** Register a field for output during write
     *
     * \param field Name of the field to be output
     */
    inline void add_output_field(const std::string field)
    {
        output_fields_.insert(field);
    }

    /** Open a database for writing time series data
     *
     *  \param output_db Pathname to the output ExodusII database
     *  \return A valid file handle for use with write_database
     *
     *  \sa write_database, write_timesteps
     */
    size_t open_database(std::string output_db);

    /** Write time series data to an open database
     *
     *  \param fh Valid file handle
     *  \param time Time to write
     *
     *  \sa open_database, write_timesteps
     */
    void write_database(size_t fh, double time);

    /** Write the Exodus results database with modifications
     *
     *  \param output_db Pathname to the output ExodusII database
     *  \param time Timestep to write
     *
     *  \sa write_database_with_fields
     */
    void write_database(std::string output_db, double time=0.0);

    /** Write database with restart fields
     *
     *  Copies the restart data fields from the input Exodus database to the
     *  output database.
     *
     *  \param output_db Pathname to the output ExodusII database
     */
    void write_database_with_fields(std::string output_db);

    /** Write time-history to database
     *
     *  This method accepts a functor that takes one integer argument (timestep)
     *  and returns the time (double) that must be written to the database. The
     *  functor should update the fields that are being written to the database.
     *  An example would be to simulate mesh motion by updating the
     *  mesh_displacement field at every timestep.
     *
     * The following example shows the use with a C++ lambda function:
     *
     *  ```
     *  double deltaT = 0.01;  // Timestep size
     *
     *  write_timesteps("inflow_history.exo", 100,
     *    [&](int tstep) {
     *        double time = tstep * deltaT;
     *
     *        // Update velocity and coordinates
     *
     *        return time;
     *  });
     *  ```
     */
    template <typename Functor, typename FieldSetType>
    void write_timesteps(
        std::string output_db,
        int num_steps,
        const FieldSetType& field_list,
        Functor lambdaFunc)
    {
        size_t fh = stkio_.create_output_mesh(output_db, stk::io::WRITE_RESULTS);
        for (auto fname: field_list) {
            stk::mesh::FieldBase* fld = stk::mesh::get_field_by_name(fname, meta_);
            if (fld != nullptr) {
                stkio_.add_field(fh, *fld, fname);
            }
        }

        for (int tstep=0; tstep < num_steps; tstep++) {
            double time = lambdaFunc(tstep);
            write_database(fh, time);
        }
    }

    template<typename Functor>
    void write_timesteps(std::string output_db, int num_steps, Functor lambdaFunc)
    {
        write_timesteps(output_db, num_steps, output_fields_, lambdaFunc);
    }

    /** Calculate the bounding box of the mesh
     *
     *  The selector can pick parts that are not contiguous. However, the
     *  bounding box returned will be the biggest box that encloses all parts
     *  selected.
     *
     *  \param selector An instance of stk::mesh::Selector to filter parts of
     *  the mesh where bounding box is calculated.
     *  \param verbose If true, then print out the bounding box to standard output.
     *
     *  \return An stk::search::Box instance containing the min and max points (3-D).
     */
    BoxType calc_bounding_box(const stk::mesh::Selector, bool verbose=true);

    /** Set automatic mesh decomposition property
     *
     *  \param decompType The decomposition type
     *
     *  Valid decomposition types are: rcb, rib, block, linear
     */
    inline void set_decomposition_type(std::string decompType)
    {
        stkio_.property_add(Ioss::Property("DECOMPOSITION_METHOD", decompType));
    }

    /** Join output/restart files when writing to disk
     *
     */
    inline void set_auto_join()
    {
        stkio_.property_add(Ioss::Property("COMPOSE_RESULTS", "YES"));
        stkio_.property_add(Ioss::Property("COMPOSE_RESTART", "YES"));
    }

    /** Force output database to use 8-bit integers
     */
    inline void set_64bit_flags()
    {
        stkio_.property_add(Ioss::Property("INTEGER_SIZE_DB",8));
        stkio_.property_add(Ioss::Property("INTEGER_SIZE_API",8));
    }

    //! Flag indicating whether the DB has been modified
    inline bool db_modified()
    {
        return db_modified_ || (output_fields_.size() > 0);
    }

    //! Force output of the results DB
    inline void set_write_flag(bool flag = true)
    {
        db_modified_ = flag;
    }

    //! Return a reference to the registered output fields
    inline const std::unordered_set<std::string>& output_fields()
    { return output_fields_; }

private:
    CFDMesh() = delete;
    CFDMesh(const CFDMesh&) = delete;
    //! Instance of the parallel MPI_COMM object
    stk::ParallelMachine& comm_;

    //! STK Metadata object
    stk::mesh::MetaData meta_;

    //! STK Bulkdata object
    stk::mesh::BulkData bulk_;

    //! Filename for the mesh input database
    std::string input_db_{""};

    //! STK Mesh I/O handler
    stk::io::StkMeshIoBroker stkio_;

    //! List of fields registered for output
    std::unordered_set<std::string> output_fields_;

    //! Flag indicating whether a results database must be written
    bool db_modified_{false};
};

}  // nalu
}  // sierra

#endif /* CFDMESH_H */
