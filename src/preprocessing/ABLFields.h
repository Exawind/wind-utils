#ifndef ABLFIELDS_H
#define ABLFIELDS_H

#include "PreProcessingTask.h"

namespace sierra {
namespace nalu {

class ABLFields: public PreProcessingTask
{
public:
    template<typename T>
    using Array2D = std::vector<std::vector<T>>;

    ABLFields(CFDMesh&, const YAML::Node&);

    virtual ~ABLFields() {}

    void initialize();

    void run();

private:
    ABLFields() = delete;
    ABLFields(const ABLFields&) = delete;

    void load(const YAML::Node&);

    void load_velocity_info(const YAML::Node&);

    void load_temperature_info(const YAML::Node&);

    void init_velocity_field();

    void init_temperature_field();

    //! STK Metadata object
    stk::mesh::MetaData& meta_;

    //! STK Bulkdata object
    stk::mesh::BulkData& bulk_;

    //! Parts of the fluid mesh where wall distance is calculated
    stk::mesh::PartVector fluid_parts_;

    std::vector<double> vHeights_;

    Array2D<double> velocity_;

    std::vector<double> THeights_;

    std::vector<double> TValues_;

    //! Dimensionality of the mesh
    int ndim_;

    bool doVelocity_;

    bool doTemperature_;
};

}
}

#endif /* ABLFIELDS_H */
