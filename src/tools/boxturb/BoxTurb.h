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

#ifndef BOXTURB_H
#define BOXTURB_H

#include <vector>
#include <memory>

#include "struct_grid/StructGrid.h"
#include "struct_grid/StructGridIx.h"
#include "struct_grid/StructField.h"

#include "core/YamlUtils.h"

namespace sierra {
namespace nalu {

class BoxTurb
{
public:
    static constexpr int ndim = StructBox::ndim;
    using Layout = sgix::LeftLayout;

    BoxTurb() = default;

    ~BoxTurb() = default;

    void load(const YAML::Node&);

    void run(const YAML::Node&);

    void write_netcdf(const std::string outfile);

    double boxVol() const
    { return boxlen_[0] * boxlen_[1] * boxlen_[2]; }

private:
    struct NCBoxTurb
    {
        int ncid;

        int sDim, xDim, yDim, zDim;
        int uid, vid, wid;
        int blenid, dxid;

        int scaleid, divCorrid;
    };

    void correct_divU(const YAML::Node&);

    void write_netcdf_dims(NCBoxTurb&);

    void write_netcdf_data(NCBoxTurb&);

#ifdef ENABLE_HYPRE
    void solve_divU_hypre(const YAML::Node&);
#endif

    void exchange_ghosts_single(BoxField<double>&);

    void exchange_ghosts(BoxField<double>&);

    void project_velocity();

    StructGrid grid_;

    std::unique_ptr<BoxField<double>> uvel_;
    std::unique_ptr<BoxField<double>> vvel_;
    std::unique_ptr<BoxField<double>> wvel_;
    std::unique_ptr<BoxField<double>> pressure_;

    std::string source_;

    double boxlen_[ndim];

    double dx_[ndim];

    double scale_factors_[ndim]{1.0, 1.0, 1.0};

    bool correct_divU_{false};

    bool apply_scaling_{false};

    friend class BoxTurbIO;
};


}  // nalu
}  // sierra


#endif /* BOXTURB_H */
