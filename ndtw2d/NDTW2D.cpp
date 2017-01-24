
#include "NDTW2D.h"

#include <Shards_BasicTopologies.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/BroadcastArg.hpp>

#include <stk_mesh/base/FindRestriction.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/Stencils.hpp>
#include <stk_mesh/base/TopologyDimensions.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>

#include <stk_io/IossBridge.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <Ionit_Initializer.h>

#include <cmath>


namespace sierra {
namespace nalu {

typedef stk::mesh::Field<double, stk::mesh::Cartesian> VectorFieldType;
typedef stk::mesh::Field<double> ScalarFieldType;

NDTW2D::NDTW2D
(
    stk::ParallelMachine& comm,
    const YAML::Node& node
):
    ndim_(2),
    comm_(comm),
    meta_(),
    bulk_(meta_,comm),
    input_db_(),
    output_db_(),
    wall_dist_name_("NDTW"),
    fluid_parts_(0),
    wall_parts_(0),
    stkIo_()
{
    load(node);
}

void NDTW2D::load(const YAML::Node& node)
{
    const YAML::Node& wdist = node["calc_ndtw2d"];

    std::cout << "Parsing input file... " << std::endl;
    input_db_ = wdist["input_db"].as<std::string>();
    output_db_ = wdist["output_db"].as<std::string>();
    auto fluid_partnames = wdist["fluid_parts"].as<std::vector<std::string>>();
    auto wall_partnames = wdist["wall_parts"].as<std::vector<std::string>>();

    if(wdist["wall_dist_name"]) {
        wall_dist_name_ = wdist["wall_dist_name"].as<std::string>();
    }

    std::cout << "Loading input mesh database... " << std::endl;
    stkIo_.reset(new stk::io::StkMeshIoBroker(comm_));
    stkIo_->add_mesh_database(input_db_, stk::io::READ_MESH);
    stkIo_->set_bulk_data(bulk_);
    stkIo_->create_input_mesh();
    stkIo_->add_all_mesh_fields_as_input_fields();

    fluid_parts_.resize(fluid_partnames.size());
    wall_parts_.resize(wall_partnames.size());

    for(int i=0; i<fluid_partnames.size(); i++) {
        stk::mesh::Part* part = meta_.get_part(fluid_partnames[i]);
        if (NULL == part) {
            throw std::runtime_error("Missing fluid part in mesh database: " +
                                     fluid_partnames[i]);
        } else {
            fluid_parts_[i] = part;
        }
    }
    for(int i=0; i<wall_partnames.size(); i++) {
        stk::mesh::Part* part = meta_.get_part(wall_partnames[i]);
        if (NULL == part) {
            throw std::runtime_error("Missing wall part in mesh database: " +
                                     wall_partnames[i]);
        } else {
            wall_parts_[i] = part;
        }
    }
}

void NDTW2D::generate_meta()
{
    std::cout << "Registering fields to metadata..." << std::endl;
    VectorFieldType& coords = meta_.declare_field<VectorFieldType>(
        stk::topology::NODE_RANK, "coordinates");
    ScalarFieldType& ndtw = meta_.declare_field<ScalarFieldType>(
        stk::topology::NODE_RANK, wall_dist_name_);

    for(auto part: fluid_parts_) {
        stk::mesh::put_field(coords, *part, ndim_);
        stk::mesh::put_field(ndtw, *part);
    }

    std::cout << "Populating bulk data..." << std::endl;
    stkIo_->populate_bulk_data();
}

void NDTW2D::calc_ndtw()
{
    stk::mesh::Selector fluid_union = stk::mesh::selectUnion(fluid_parts_);
    stk::mesh::Selector wall_union = stk::mesh::selectUnion(wall_parts_);

    const stk::mesh::BucketVector& fluid_bkts = bulk_.get_buckets(
        stk::topology::NODE_RANK, fluid_union);
    const stk::mesh::BucketVector& wall_bkts = bulk_.get_buckets(
        stk::topology::NODE_RANK, wall_union);

    VectorFieldType* coords = meta_.get_field<VectorFieldType>(
        stk::topology::NODE_RANK, "coordinates");
    ScalarFieldType* ndtw = meta_.get_field<ScalarFieldType>(
        stk::topology::NODE_RANK, wall_dist_name_);

    for(size_t ib=0; ib < fluid_bkts.size(); ib++) {
        stk::mesh::Bucket& fbkt = *fluid_bkts[ib];
        double* xyz = stk::mesh::field_data(*coords, fbkt);
        double* wdist = stk::mesh::field_data(*ndtw, fbkt);

        for(size_t in=0; in < fbkt.size(); in++) {
            double min_dist = std::numeric_limits<double>::max();
            for(size_t jb=0; jb < wall_bkts.size(); jb++) {
                stk::mesh::Bucket& wbkt = *wall_bkts[jb];
                double* wxyz = stk::mesh::field_data(*coords, wbkt);

                for(size_t jn=0; jn< wbkt.size(); jn++) {
                    double dist_calc = 0.0;
                    for(int j=0; j<ndim_; j++) {
                        double dst = xyz[in*ndim_+j] - wxyz[jn*ndim_+j];
                        dist_calc += dst * dst;
                    }
                    if (dist_calc < min_dist) min_dist = dist_calc;
                }
            }
            wdist[in] = std::sqrt(min_dist);
        }
    }
}

void NDTW2D::run()
{
    generate_meta();

    size_t fh = stkIo_->create_output_mesh(output_db_, stk::io::WRITE_RESULTS);
    ScalarFieldType* ndtw = meta_.get_field<ScalarFieldType>(
        stk::topology::NODE_RANK, wall_dist_name_);
    stkIo_->add_field(fh, *ndtw);

    // Perform actual calculation
    calc_ndtw();

    // Dump the wall distance into the output mesh database
    stkIo_->begin_output_step(fh, 0.0);
    stkIo_->write_defined_output_fields(fh);
    stkIo_->end_output_step(fh);
}


}  // nalu
}  // sierra
