#ifndef NDTW2D_H
#define NDTW2D_H

#include "PreProcessingTask.h"

namespace sierra {
namespace nalu {

/** 2-D Nearest distance to wall calculator
 *
 * Calculates a new field NDTW containing the wall distance for 2-D airfoil-like
 * applications used in RANS wall models.
 */
class NDTW2D: public PreProcessingTask
{
public:
    NDTW2D(CFDMesh&, const YAML::Node&);

    virtual ~NDTW2D() {}

    //! Initialize the NDTW field and register for output
    virtual void initialize();

    //! Calculate wall distance and update NDTW field
    virtual void run();

private:
    NDTW2D() = delete;
    NDTW2D(const NDTW2D&) = delete;

    //! Helper method to parse YAML file and initialize parameters
    void load(const YAML::Node&);

    //! Calculate minimum distance to wall for all nodes in the mesh
    /**
     * Currently performs brute-force search and is not suitable for 3-D
     * applications.
     *
     * \todo Move to percept when publicly available
     */
    void calc_ndtw();

    //! STK Metadata object
    stk::mesh::MetaData& meta_;

    //! STK Bulkdata object
    stk::mesh::BulkData& bulk_;

    //! Parts of the fluid mesh where wall distance is calculated
    stk::mesh::PartVector fluid_parts_;

    //! Part names of the wall boundaries
    stk::mesh::PartVector wall_parts_;

    //! Field name for wall distance (default: NDTW)
    std::string wall_dist_name_;

    //! Dimensionality of the mesh
    int ndim_;

};

} // nalu
} // sierra

#endif /* NDTW2D_H */
