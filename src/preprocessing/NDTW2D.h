#ifndef NDTW2D_H
#define NDTW2D_H

#include "PreProcessingTask.h"

namespace sierra {
namespace nalu {

class NDTW2D: public PreProcessingTask
{
public:
    NDTW2D(CFDMesh&, const YAML::Node&);

    virtual ~NDTW2D() {}

    virtual void initialize();

    virtual void run();

private:
    NDTW2D() = delete;
    NDTW2D(const NDTW2D&) = delete;

    void load(const YAML::Node&);

    void calc_ndtw();

    //! STK Metadata object
    stk::mesh::MetaData& meta_;

    //! STK Bulkdata object
    stk::mesh::BulkData& bulk_;

    //! Parts of the fluid mesh where wall distance is calculated
    stk::mesh::PartVector fluid_parts_;

    //! Part names of the wall boundaries
    stk::mesh::PartVector wall_parts_;

    //! Field name for wall distance
    std::string wall_dist_name_;

    //! Dimensionality of the mesh
    int ndim_;

};

} // nalu
} // sierra

#endif /* NDTW2D_H */
