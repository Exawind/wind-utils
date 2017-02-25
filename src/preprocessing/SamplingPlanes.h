#ifndef SAMPLINGPLANES_H
#define SAMPLINGPLANES_H

#include "PreProcessingTask.h"

namespace sierra {
namespace nalu {

class SamplingPlanes: public PreProcessingTask
{
public:
    SamplingPlanes(CFDMesh&, const YAML::Node&);

    virtual ~SamplingPlanes() {}

    virtual void initialize();

    virtual void run();

private:
    SamplingPlanes() = delete;
    SamplingPlanes(const SamplingPlanes&) = delete;

    void load(const YAML::Node&);

    //! Use fluid Realm mesh to estimate the x-y bounding box for the sampling
    //! planes.
    void calc_bounding_box();

    //! Generate entities and update coordinates for a given sampling plane
    void generate_zplane(const double);

    //! STK Metadata object
    stk::mesh::MetaData& meta_;

    //! STK Bulkdata object
    stk::mesh::BulkData& bulk_;

    //! Heights where the averaging planes are generated
    std::vector<double> heights_;

    //! Bounding box of the original mesh
    std::array<std::array<double,3>,2> bBox_;

    //! Format specification for the part name
    std::string name_format_;

    //! Fluid realm part (to determine mesh bounding box)
    std::string fluidPart_;

    //! Spatial resolution in x and y directions
    double dx_;
    double dy_;

    //! Number of nodes in x and y directions
    size_t nx_, ny_;

    //! Number of elements in x and y directions
    size_t mx_, my_;

    //! Dimensionality of the mesh
    int ndim_;
};

}  // nalu
}  // sierra

#endif /* SAMPLINGPLANES_H */
