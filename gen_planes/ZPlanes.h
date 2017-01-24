#ifndef ZPLANES_H
#define ZPLANES_H

#include "stk_util/parallel/Parallel.hpp"

#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/Entity.hpp"
#include "stk_mesh/base/Field.hpp"

#include "stk_io/IossBridge.hpp"
#include "stk_io/StkMeshIoBroker.hpp"

#include "yaml-cpp/yaml.h"

#include <vector>
#include <memory>
#include <array>

namespace sierra {
namespace nalu {

class ZPlanes
{
public:
    ZPlanes(
        stk::ParallelMachine&,
        const YAML::Node&);

    /**
     * Generate sampling planes and output new mesh database
     */
    void run();

private:
    ZPlanes();
    ZPlanes(const ZPlanes&);

    //! Load user inputs from YAML input file
    void load(const YAML::Node&);

    //! Use fluid Realm mesh to estimate the x-y bounding box for the sampling
    //! planes.
    void calc_bounding_box();

    //! Register new parts with mesh MetaData for the sampling planes
    void generate_meta();

    //! Generate entities and update coordinates for a given sampling plane
    void generate_zplane(const double);

    //! Spatial dimension (always 3)
    const int ndim_;

    //! Instance of the parallel MPI_COMM for use with StkMeshIoBroker
    stk::ParallelMachine& comm_;

    //! STK Metadata object
    stk::mesh::MetaData meta_;

    //! STK Bulkdata object
    stk::mesh::BulkData bulk_;

    //! Filename of the input mesh (must contain the Fluid Realm mesh)
    std::string input_mesh_;

    //! Filename of the output mesh.
    //!
    //! Writes the contents of the input mesh as well as the new planes
    std::string output_mesh_;

    //! Heights where the averaging planes are generated
    std::vector<double> heights_;

    //! Format specification for the part name
    std::string name_format_;

    //! Spatial resolution in x and y directions
    double dx_;
    double dy_;

    //! Number of nodes in x and y directions
    size_t nx_, ny_;

    //! Number of elements in x and y directions
    size_t mx_, my_;


    //! STK Mesh I/O handler
    std::unique_ptr<stk::io::StkMeshIoBroker> stkIo_;

    //! Fluid realm part (to determine mesh bounding box)
    std::string fluidPart_;

    //! Bounding box of the original mesh
    std::array<std::array<double,3>,2> bBox_;

};

} // namespace nalu
} // namespace sierra

#endif /* ZPLANES_H */
