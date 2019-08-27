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

#include <cassert>
#include <cmath>

#include "tools/boxturb/BoxTurb.h"
#include "tools/boxturb/BoxTurbIO.h"
#include "core/ParallelInfo.h"

#include "struct_grid/StructGridUtils.h"
#include "struct_grid/FVStencil.h"
#include "struct_grid/NeighborMap.h"

#ifdef ENABLE_HYPRE
#include "struct_grid/HypreStructSolver.h"
#endif

#include "netcdf.h"
#include "netcdf_par.h"

namespace sierra {
namespace nalu {

namespace {

inline void
check_nc_error(int ierr)
{
    if (ierr != NC_NOERR)
        throw std::runtime_error(
            "NetCDF Error: " + std::string(nc_strerror(ierr)));
}

#define BTNC_CALL(funcarg)                                                     \
    do {                                                                       \
        int ierr = funcarg;                                                    \
        check_nc_error(ierr);                                                  \
    } while (0)

}

void BoxTurb::load(const YAML::Node& node)
{
    using idx_t = SGTraits::idx_t;
    const auto& pinfo = get_mpi();

    wind_utils::get_optional(node, "correct_divergence", correct_divU_);

    auto box_dims = node["box_dims"].as<std::vector<int>>();
    auto box_len = node["box_len"].as<std::vector<double>>();
    assert(box_dims.size() == SGTraits::ndim);
    assert(box_len.size() == SGTraits::ndim);

    if (correct_divU_)
        grid_.set_num_ghost(1);

    grid_.set_global_grid(box_dims[0], box_dims[1], box_dims[2]);
    if (node["partitions"]) {
        auto parts = node["partitions"].as<std::vector<int>>();
        grid_.set_partitions(parts[0], parts[1], parts[2]);
    } else
        grid_.set_partitions(pinfo.size());

    const auto* gsize = grid_.global().size;
    for (int i = 0; i < SGTraits::ndim; ++i) {
        boxlen_[i] = box_len[i];
        dx_[i] = box_len[i] / static_cast<double>(gsize[i]);
    }

    const auto& local = grid_.local();
    uvel_.reset(new BoxField<double>(local));
    vvel_.reset(new BoxField<double>(local));
    wvel_.reset(new BoxField<double>(local));

    const std::string ftype("windsim");
    BoxTurbIO turbIO(*this);
    turbIO.load(ftype, node);

    // Scaling factors
    std::string scale_type = "none";
    wind_utils::get_optional(node, "scale_type", scale_type);

    if (scale_type != "none") {
        const auto sfac = node["scaling_factors"].as<std::vector<double>>();
        std::copy(sfac.begin(), sfac.end(), std::begin(scale_factors_));

        // Assume user wants to apply scaling
        apply_scaling_ = true;
    }

    wind_utils::get_optional(node, "apply_scaling", apply_scaling_);
}

void BoxTurb::run(const YAML::Node& node)
{
    using idx_t = SGTraits::idx_t;
    if (apply_scaling_) {
        sgix::ijk_loop(
            grid_.local(),
            [&](idx_t i, idx_t j, idx_t k) {
                (*uvel_)(i, j, k) *= scale_factors_[0];
                (*vvel_)(i, j, k) *= scale_factors_[1];
                (*wvel_)(i, j, k) *= scale_factors_[2];
            });
    }

    if (correct_divU_) correct_divU(node);

    if (node["output"]) {
        const auto output = node["output"].as<std::string>();
        write_netcdf(output);
    }
}

void BoxTurb::write_netcdf(const std::string outfile)
{
    const auto& pinfo = get_mpi();
    int ierr = 0;
    BoxTurb::NCBoxTurb bt;

    pinfo.info() << "Begin output in NetCDF format: " << outfile << std::endl;

    ierr = nc_create_par(
        outfile.c_str(), NC_CLOBBER | NC_NETCDF4 | NC_MPIIO,
        pinfo.comm(), MPI_INFO_NULL, &bt.ncid);
    check_nc_error(ierr);

    write_netcdf_dims(bt);
    write_netcdf_data(bt);

    ierr = nc_close(bt.ncid);
    check_nc_error(ierr);

    pinfo.info() << "NetCDF file written successfully: " << outfile << std::endl;
}

void BoxTurb::write_netcdf_dims(BoxTurb::NCBoxTurb& bt)
{
    int ierr = 0;

    const auto* nx = grid_.global().size;

    // Define dimensions for NetCDF file
    ierr = nc_def_dim(bt.ncid, "ndim", static_cast<int>(ndim), &bt.sDim);
    ierr = nc_def_dim(bt.ncid, "nx", nx[0], &bt.xDim);
    ierr = nc_def_dim(bt.ncid, "ny", nx[1], &bt.yDim);
    ierr = nc_def_dim(bt.ncid, "nz", nx[2], &bt.zDim);

    // Flags
    ierr = nc_def_var(bt.ncid, "divergence_correction", NC_INT, 0, NULL, &bt.divCorrid);

    // Grid resolution information
    ierr = nc_def_var(bt.ncid, "box_lengths", NC_DOUBLE, 1, &bt.sDim, &bt.blenid);
    ierr = nc_def_var(bt.ncid, "dx", NC_DOUBLE, 1, &bt.sDim, &bt.dxid);

    // Velocity scaling factors
    ierr = nc_def_var(bt.ncid, "scaling_factors", NC_DOUBLE, 1, &bt.sDim, &bt.scaleid);

    // Define velocity arrays
    const std::vector<int> threeD{bt.xDim, bt.yDim, bt.zDim};
    ierr = nc_def_var(
        bt.ncid, "uvel", NC_DOUBLE, static_cast<int>(ndim), threeD.data(), &bt.uid);
    ierr = nc_def_var(
        bt.ncid, "vvel", NC_DOUBLE, static_cast<int>(ndim), threeD.data(), &bt.vid);
    ierr = nc_def_var(
        bt.ncid, "wvel", NC_DOUBLE, static_cast<int>(ndim), threeD.data(), &bt.wid);

    // Attributes
    const std::string vUnit = "m/s";
    const std::string lUnit = "m";
    ierr = nc_put_att_text(bt.ncid, NC_GLOBAL, "title", source_.size()+1, source_.c_str());

    ierr = nc_put_att_text(
        bt.ncid, bt.dxid, "units", lUnit.size() + 1, lUnit.c_str());
    ierr = nc_put_att_text(
        bt.ncid, bt.blenid, "units", lUnit.size() + 1, lUnit.c_str());
    ierr = nc_put_att_text(
        bt.ncid, bt.uid, "units", vUnit.size() + 1, vUnit.c_str());
    ierr = nc_put_att_text(
        bt.ncid, bt.vid, "units", vUnit.size() + 1, vUnit.c_str());
    ierr = nc_put_att_text(
        bt.ncid, bt.wid, "units", vUnit.size() + 1, vUnit.c_str());

    ierr = nc_enddef(bt.ncid);
    check_nc_error(ierr);
}

void BoxTurb::write_netcdf_data(BoxTurb::NCBoxTurb& bt)
{
    using idx_t = SGTraits::idx_t;
    const auto& pinfo = get_mpi();
    int ierr = 0;

    if (pinfo.master()) {
        nc_put_var_double(bt.ncid, bt.blenid, boxlen_);
        nc_put_var_double(bt.ncid, bt.dxid, dx_);
        nc_put_var_double(bt.ncid, bt.scaleid, scale_factors_);

        int divcorr = (correct_divU_) ? 1 : 0;
        nc_put_var(bt.ncid, bt.divCorrid, &divcorr);
    }

    const auto rbox = sgix::real_box(grid_.local());
    BoxField<double> buffer(rbox);
    std::vector<size_t> start(SGTraits::ndim), count(SGTraits::ndim);
    start[0] = rbox.start[0];
    start[1] = rbox.start[1];
    start[2] = rbox.start[2];
    count[0] = rbox.size[0];
    count[1] = rbox.size[1];
    count[2] = rbox.size[2];

    {
        const auto& uvel = *uvel_;
        sgix::ijk_loop(rbox, [&](idx_t i, idx_t j, idx_t k) {
            buffer(i, j, k) = uvel(i, j, k);
        });

        ierr = nc_put_vara_double(
            bt.ncid, bt.uid, start.data(), count.data(), buffer.data());
        check_nc_error(ierr);
    }

    {
        const auto& vvel = *vvel_;
        sgix::ijk_loop(rbox, [&](idx_t i, idx_t j, idx_t k) {
            buffer(i, j, k) = vvel(i, j, k);
        });

        ierr = nc_put_vara_double(
            bt.ncid, bt.vid, start.data(), count.data(), buffer.data());
        check_nc_error(ierr);
    }

    {
        const auto& wvel = *wvel_;
        sgix::ijk_loop(rbox, [&](idx_t i, idx_t j, idx_t k) {
            buffer(i, j, k) = wvel(i, j, k);
        });

        ierr = nc_put_vara_double(
            bt.ncid, bt.wid, start.data(), count.data(), buffer.data());
        check_nc_error(ierr);
    }
}

void BoxTurb::correct_divU(const YAML::Node& node)
{
#ifdef ENABLE_HYPRE
    pressure_.reset(new BoxField<double>(grid_.local()));

    // Only HYPRE solver for now
    solve_divU_hypre(node);

    exchange_ghosts(*pressure_);
    project_velocity();
#else
    throw std::runtime_error("Divergence correction requires HYPRE support");
#endif
}

#ifdef ENABLE_HYPRE
void BoxTurb::solve_divU_hypre(const YAML::Node& node)
{
    using idx_t = SGTraits::idx_t;
    HypreStructSolver solver(grid_);
    const auto& snode = node["solver_settings"];

    solver.init_hypre_grid();
    solver.init_solver(snode);

    std::vector<double> stencil_entries(7);
    const double dVol = dx_[0] * dx_[1] * dx_[2];
    const double dx2 = 1.0 / (dx_[0] * dx_[0]);
    const double dy2 = 1.0 / (dx_[1] * dx_[1]);
    const double dz2 = 1.0 / (dx_[2] * dx_[2]);
    const double dxdy = dx_[0] * dx_[1];
    const double dydz = dx_[1] * dx_[2];
    const double dxdz = dx_[0] * dx_[2];

    stencil_entries[fvm::CENTER] =  2.0 * dVol * ( dx2 + dy2 + dz2);
    stencil_entries[fvm::WEST]   = -dVol * dx2;
    stencil_entries[fvm::EAST]   = -dVol * dx2;
    stencil_entries[fvm::SOUTH]  = -dVol * dy2;
    stencil_entries[fvm::NORTH]  = -dVol * dy2;
    stencil_entries[fvm::BOTTOM] = -dVol * dz2;
    stencil_entries[fvm::TOP]    = -dVol * dz2;

    solver.populate_matrix(stencil_entries);

    const auto rbox = sgix::real_box(grid_.local());
    const auto ncells = sgix::num_cells(rbox);
    std::vector<double> values(ncells);
    const auto& uvel = *uvel_;
    const auto& vvel = *vvel_;
    const auto& wvel = *wvel_;

    SGTraits::size_t m = 0;
    sgix::kji_loop(
        rbox,
        [&](idx_t i, idx_t j, idx_t k) {
            values[m] = - 0.5 * (
                (uvel(i-1, j, k) - uvel(i+1, j, k)) * dydz +
                (vvel(i, j+1, k) - vvel(i, j-1, k)) * dxdz +
                (wvel(i, j, k+1) - wvel(i, j, k-1)) * dxdy);
            m++;
        });

    solver.populate_rhs(values);
    solver.solve();
    solver.get_solution(values);

    m = 0;
    auto& pres = *pressure_;
    sgix::kji_loop(
        rbox,
        [&](idx_t i, idx_t j, idx_t k) {
            pres(i, j, k) = values[m++];
        });
}
#endif

void BoxTurb::exchange_ghosts_single(BoxField<double>& field)
{
    using idx_t = SGTraits::idx_t;
    const auto partitions = grid_.partitions();
    const auto& global = grid_.global();
    const auto& rbox = sgix::real_box(grid_.local());

    if (partitions[2] == 1) {
        const idx_t kmin = global.start[2];
        const idx_t kmax = global.end[2] - 1;

        for(idx_t i=rbox.start[0]; i < rbox.end[0]; ++i)
            for (idx_t j=rbox.start[1]; j < rbox.end[1]; ++j) {
                field(i, j, kmin - 1) = field(i, j, kmax);
                field(i, j, kmax + 1) = field(i, j, kmin);
            }
    }

    if (partitions[1] == 1) {
        const idx_t jmin = global.start[1];
        const idx_t jmax = global.end[1] - 1;

        for (idx_t i = rbox.start[0]; i < rbox.end[0]; ++i)
            for (idx_t k = rbox.start[2]; k < rbox.end[2]; ++k) {
                field(i, jmin - 1, k) = field(i, jmax, k);
                field(i, jmax + 1, k) = field(i, jmin, k);
            }
    }

    if (partitions[0] == 1) {
        const idx_t imin = global.start[0];
        const idx_t imax = global.end[0] - 1;

        for (idx_t j = rbox.start[1]; j < rbox.end[1]; ++j)
            for (idx_t k = rbox.start[2]; k < rbox.end[2]; ++k) {
                field(imin - 1, j, k) = field(imax, j, k);
                field(imax + 1, j, k) = field(imin, j, k);
            }
    }
}

void BoxTurb::exchange_ghosts(BoxField<double>& field)
{
    static constexpr int ndim = SGTraits::ndim;
    static constexpr int nfaces = 2 * ndim;
    using idx_t = SGTraits::idx_t;

    const int myid = get_mpi().rank();
    const NeighborMap<sgix::RightLayout> mpiMap(grid_);
    std::vector<StructBox> ghost_box(nfaces);
    std::vector<std::vector<double>> sendbuf(nfaces), recvbuf(nfaces);

    const auto local = grid_.local();
    const auto* data = field.data();
    const auto idxp = sgix::PeriodicIndexer<sgix::LeftLayout>(local);

    idx_t ioff[ndim] = {0, 0, 0};
    // Identify ghost layers and populate send buffers
    for (int f=0; f < nfaces; ++f) {
        const int id = f / 2;
        ghost_box[f] = sgix::ghost_layer(
            local, static_cast<fvm::FVStencil>(f + 1));

        const auto ncells = sgix::num_cells(ghost_box[f]);
        sendbuf[f].resize(ncells);
        recvbuf[f].resize(ncells);

        for (int i =0; i < ndim; i++) ioff[i] = 0;
        ioff[id] = ((f % 2) == 0)
            ? -local.nghost[id] - 1
            :  local.nghost[id] + 1;

        SGTraits::size_t m = 0;
        sgix::ijk_loop(
            ghost_box[f],
            [&](idx_t i, idx_t j, idx_t k) {
                const auto ii = idxp(i + ioff[0], j + ioff[1], k + ioff[2]);
                sendbuf[f][m] = data[ii];
                m++;
            });
    }

    // Send/recv
    const auto comm = get_mpi().comm();
    std::vector<MPI_Request> request(2*nfaces);
    std::vector<MPI_Status> status(2*nfaces);
    for (int f=0; f < nfaces; ++f) {
        const int neighbor = mpiMap(static_cast<fvm::FVStencil>(f+1));
        const int tag = 100 + f;

        MPI_Irecv(
            recvbuf[f].data(), recvbuf[f].size(), MPI_DOUBLE, neighbor, tag,
            comm, &request[2 * f]);

        MPI_Isend(
            sendbuf[f].data(), sendbuf[f].size(), MPI_DOUBLE, neighbor, tag,
            comm, &request[2 * f + 1]);
    }

    MPI_Waitall(2*nfaces, request.data(), status.data());

    for (int f=0; f < nfaces; ++f) {
        SGTraits::size_t m = 0;
        sgix::ijk_loop(
            ghost_box[f],
            [&](idx_t i, idx_t j, idx_t k) {
                field(i, j, k) = recvbuf[f][m];
                m++;
        });
    }
}

void BoxTurb::project_velocity()
{
    using idx_t = SGTraits::idx_t;
    auto& uvel = *uvel_;
    auto& vvel = *vvel_;
    auto& wvel = *wvel_;
    auto& pres = *pressure_;

    const auto inv_dx = 0.5 / dx_[0];
    const auto inv_dy = 0.5 / dx_[1];
    const auto inv_dz = 0.5 / dx_[2];

    const auto rbox = sgix::real_box(grid_.local());
    sgix::ijk_loop(
        rbox,
        [&](idx_t i, idx_t j, idx_t k) {
            uvel(i, j, k) += (pres(i+1, j, k) - pres(i-1, j, k)) * inv_dx;
            vvel(i, j, k) -= (pres(i, j+1, k) - pres(i, j-1, k)) * inv_dy;
            wvel(i, j, k) -= (pres(i, j, k+1) - pres(i, j, k-1)) * inv_dz;
        });

#if 0
    exchange_ghosts(*uvel_);
    exchange_ghosts(*vvel_);
    exchange_ghosts(*wvel_);

    double res = 0.0;
    const double dVol = dx_[0] * dx_[1] * dx_[2];
    const double dxdy = dx_[0] * dx_[1];
    const double dydz = dx_[1] * dx_[2];
    const double dxdz = dx_[0] * dx_[2];
    sgix::ijk_loop(
        rbox,
        [&](idx_t i, idx_t j, idx_t k) {
            double tmp = (uvel(i-1, j, k) - uvel(i+1, j, k)) * dydz +
                (vvel(i, j+1, k) - vvel(i, j-1, k)) * dxdz +
                (wvel(i, j, k+1) - wvel(i, j, k-1)) * dxdy;
            res += tmp * 0.5;
    });

    double gres = 0.0;
    MPI_Reduce(&res, &gres, 1, MPI_DOUBLE, MPI_SUM, 0, get_mpi().comm());
    std::cout << gres / boxVol() << std::endl;
#endif
}

}  // nalu
}  // sierra
