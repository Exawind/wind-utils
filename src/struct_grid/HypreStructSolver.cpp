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

#include <cctype>
#include <algorithm>

#include "struct_grid/HypreStructSolver.h"
#include "struct_grid/StructGridIx.h"
#include "struct_grid/FVStencil.h"
#include "core/ParallelInfo.h"

namespace sierra {
namespace nalu {

HypreStructSolver::~HypreStructSolver()
{
    if (systemInitialized_) {
        HYPRE_StructGridDestroy(hgrid_);
        HYPRE_StructStencilDestroy(stencil_);
        HYPRE_StructMatrixDestroy(Amat_);
        HYPRE_StructVectorDestroy(rhs_);
        HYPRE_StructVectorDestroy(sln_);
        systemInitialized_ = false;
    }

    if (solverInitialized_)
        solverDestroyPtr_(solver_);
}

void HypreStructSolver::init_solver(const YAML::Node& node)
{
    std::string method = "SMG";
    wind_utils::get_optional(node, "method", method);
    std::transform(method.begin(), method.end(), method.begin(),
                   [](unsigned char c) { return std::toupper(c); });
    method_ = method;

    std::string precond = "none";
    wind_utils::get_optional(node, "preconditioner", precond);

    wind_utils::get_optional(node, "max_iterations", maxIterations_);
    wind_utils::get_optional(node, "tolerance", tolerance_);
    wind_utils::get_optional(node, "print_level", printLevel_);
    wind_utils::get_optional(node, "log_level", logLevel_);

    if (method == "SMG")
        setup_smg_solver(node);
    else if (method == "PFMG")
        setup_pfmg_solver(node);
    else
        throw std::runtime_error("Invalid HYPRE solver specified: " + method);
}

void HypreStructSolver::setup_smg_solver(const YAML::Node&)
{
    const auto comm = get_mpi().comm();

    HYPRE_StructSMGCreate(comm, &solver_);

    HYPRE_StructSMGSetMaxIter(solver_, maxIterations_);
    HYPRE_StructSMGSetTol(solver_, tolerance_);
    HYPRE_StructSMGSetPrintLevel(solver_, printLevel_);
    HYPRE_StructSMGSetLogging(solver_, logLevel_);

    solverDestroyPtr_ = &HYPRE_StructSMGDestroy;
    solverSetupPtr_ = &HYPRE_StructSMGSetup;
    solverSolvePtr_ = &HYPRE_StructSMGSolve;
    solverPrecondPtr_ = nullptr;
    solverItersPtr_ = &HYPRE_StructSMGGetNumIterations;
    solverResNormPtr_ = &HYPRE_StructSMGGetFinalRelativeResidualNorm;

    usePrecond_ = false;
    solverInitialized_ = true;
}

void HypreStructSolver::setup_pfmg_solver(const YAML::Node&)
{
    const auto comm = get_mpi().comm();

    HYPRE_StructPFMGCreate(comm, &solver_);

    HYPRE_StructPFMGSetMaxIter(solver_, maxIterations_);
    HYPRE_StructPFMGSetTol(solver_, tolerance_);
    HYPRE_StructPFMGSetPrintLevel(solver_, printLevel_);
    HYPRE_StructPFMGSetLogging(solver_, logLevel_);

    solverDestroyPtr_ = &HYPRE_StructPFMGDestroy;
    solverSetupPtr_ = &HYPRE_StructPFMGSetup;
    solverSolvePtr_ = &HYPRE_StructPFMGSolve;
    solverPrecondPtr_ = nullptr;
    solverItersPtr_ = &HYPRE_StructPFMGGetNumIterations;
    solverResNormPtr_ = &HYPRE_StructPFMGGetFinalRelativeResidualNorm;

    usePrecond_ = false;
    solverInitialized_ = true;
}

void HypreStructSolver::init_hypre_grid()
{
    const auto& pinfo = get_mpi();

    // Create hypre structured grid definition
    HYPRE_StructGridCreate(pinfo.comm(), ndim, &hgrid_);

    // Assume only 1 box per processor for now
    const auto rbox = sgix::real_box(grid_.local());
    for (int d=0; d < ndim; d++) {
        ilower_[d] = rbox.start[d];
        iupper_[d] = rbox.end[d] - 1;
    }
    HYPRE_StructGridSetExtents(hgrid_, ilower_, iupper_);

    if (isPeriodic_) {
        const auto& global = grid_.global();
        for (int d=0; d < ndim; d++)
            periodic_[d] = global.size[d];

        HYPRE_StructGridSetPeriodic(hgrid_, periodic_);
    }

    HYPRE_StructGridSetNumGhost(hgrid_, num_ghost_);
    HYPRE_StructGridAssemble(hgrid_);

    // Create finite-volume stencil
    HYPRE_StructStencilCreate(ndim, fvm::NUM_STENCIL, &stencil_);
    HYPRE_Int offsets[ndim];
    for (int s = 0; s < fvm::NUM_STENCIL; ++s) {
        for (int d=0; d < ndim; d++)
            offsets[d] = fvm::fv_offsets[s][d];
        HYPRE_StructStencilSetElement(stencil_, s, offsets);
    }

    // Initialize matrix and rhs/solution vectors
    HYPRE_StructMatrixCreate(pinfo.comm(), hgrid_, stencil_, &Amat_);
    HYPRE_StructMatrixInitialize(Amat_);

    HYPRE_StructVectorCreate(pinfo.comm(), hgrid_, &rhs_);
    HYPRE_StructVectorInitialize(rhs_);

    HYPRE_StructVectorCreate(pinfo.comm(), hgrid_, &sln_);
    HYPRE_StructVectorInitialize(sln_);
    // Call assemble right away as we will not set this
    HYPRE_StructVectorAssemble(sln_);

    systemInitialized_ = true;
}

void HypreStructSolver::solve()
{
    const auto& pinfo = get_mpi();

    if (usePrecond_) {
        solverPrecondPtr_(
            (HYPRE_Solver) solver_, (HYPRE_PtrToSolverFcn)precondSolvePtr_,
            (HYPRE_PtrToSolverFcn)precondSetupPtr_, (HYPRE_Solver) precond_);
    }

    solverSetupPtr_(solver_, Amat_, rhs_, sln_);
    solverSolvePtr_(solver_, Amat_, rhs_, sln_);

    solverItersPtr_(solver_, &numIterations_);
    solverResNormPtr_(solver_, &finalResNorm_);
    pinfo.info() << method_
                 << ": Iters = " << numIterations_
                 << "; Rel. res. norm = " << finalResNorm_ << std::endl;
}

void HypreStructSolver::populate_matrix(const std::vector<double>& entries)
{
    using idx_t = SGTraits::idx_t;
    const int num_entries = entries.size();
    assert(num_entries == NUM_STENCIL);
    std::vector<HYPRE_Int> stencil_entries(num_entries);

    for (int d=0; d < num_entries; ++d)
        stencil_entries[d] = d;

    const auto rbox = sgix::real_box(grid_.local());
    const int ncells = sgix::num_cells(rbox);
    std::vector<double> values(ncells * num_entries);

    size_t m = 0;
    sgix::kji_loop(
        rbox,
        [&](idx_t, idx_t, idx_t) {
            for (int d=0; d < num_entries; d++) {
                values[m] = entries[d];
                m++;
            }
        });

    HYPRE_StructMatrixSetBoxValues(
        Amat_, ilower_, iupper_, num_entries, stencil_entries.data(),
        values.data());
    HYPRE_StructMatrixAssemble(Amat_);
}

void HypreStructSolver::populate_rhs(std::vector<double>& rhs)
{
    HYPRE_StructVectorSetBoxValues(rhs_, ilower_, iupper_, rhs.data());
    HYPRE_StructVectorAssemble(rhs_);
}

void HypreStructSolver::get_solution(std::vector<double>& sln)
{
    HYPRE_StructVectorGetBoxValues(sln_, ilower_, iupper_, sln.data());
}

}  // nalu
}  // sierra
