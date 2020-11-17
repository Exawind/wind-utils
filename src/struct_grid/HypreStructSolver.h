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

#ifndef HYPRESTRUCTSOLVER_H
#define HYPRESTRUCTSOLVER_H

#include <cassert>
#include "struct_grid/StructGrid.h"
#include "core/YamlUtils.h"

#include "HYPRE_struct_ls.h"
#include "HYPRE_krylov.h"

namespace sierra {
namespace nalu {

/** Interface to HYPRE Struct grid and solvers
 */
class HypreStructSolver
{
public:
    static constexpr int ndim = SGTraits::ndim;

    /**
     *  @param grid Global mesh information
     */
    HypreStructSolver(StructGrid grid) : grid_(grid) {}

    virtual ~HypreStructSolver();

    virtual void init_solver(const YAML::Node&);

    virtual void init_hypre_grid();

    virtual void populate_matrix(const std::vector<double>&);

    virtual void populate_rhs(std::vector<double>&);

    virtual void get_solution(std::vector<double>&);

    virtual void solve();

protected:
    void setup_smg_solver(const YAML::Node&);

    void setup_pfmg_solver(const YAML::Node&);

    //! The actual grid instance
    StructGrid grid_;

    //! Local indices of the lower corner of the bounding box of the block that
    //! is on the current MPI rank
    HYPRE_Int ilower_[ndim];

    //! Local indices of the upper right corner of the bounding box of the block
    //! that is on the current MPI rank
    HYPRE_Int iupper_[ndim];

    //! Index of the periodic boundary indicating periodicity in each direction
    HYPRE_Int periodic_[ndim];

    //! Number of ghost cell layers
    HYPRE_Int num_ghost_[2*ndim]{1, 1, 1, 1, 1, 1};

    //! HYPRE object representing the global mesh
    HYPRE_StructGrid hgrid_;

    //! HYPRE object representing the stencil at each cell (e.g., finite volume)
    HYPRE_StructStencil stencil_;

    //! HYPRE Struct Matrix instance
    HYPRE_StructMatrix Amat_;

    //! HYPRE RHS vector instance
    HYPRE_StructVector rhs_;

    //! HYPRE Solution vector instance
    HYPRE_StructVector sln_;

    //! HYPRE structured grid solver
    HYPRE_StructSolver solver_;

    //! HYPRE structured grid preconditioner for use with Krylov solvers
    HYPRE_StructSolver precond_;

    HYPRE_Int (*solverDestroyPtr_)(HYPRE_StructSolver);
    HYPRE_Int (*solverSetupPtr_)(
        HYPRE_StructSolver,
        HYPRE_StructMatrix,
        HYPRE_StructVector,
        HYPRE_StructVector);
    HYPRE_Int (*solverSolvePtr_)(
        HYPRE_StructSolver,
        HYPRE_StructMatrix,
        HYPRE_StructVector,
        HYPRE_StructVector);
    HYPRE_Int (*solverPrecondPtr_)(
        HYPRE_Solver,
        HYPRE_PtrToSolverFcn,
        HYPRE_PtrToSolverFcn,
        HYPRE_Solver);

    HYPRE_Int (*solverItersPtr_)(
        HYPRE_StructSolver, HYPRE_Int*);
    HYPRE_Int (*solverResNormPtr_)(
        HYPRE_StructSolver, double*);

    HYPRE_Int (*precondDestroyPtr_)(HYPRE_StructSolver);
    HYPRE_Int (*precondSetupPtr_)(
        HYPRE_StructSolver,
        HYPRE_StructMatrix,
        HYPRE_StructVector,
        HYPRE_StructVector);
    HYPRE_Int (*precondSolvePtr_)(
        HYPRE_StructSolver,
        HYPRE_StructMatrix,
        HYPRE_StructVector,
        HYPRE_StructVector);

    //! Solver used for solving the linear system
    std::string method_;

    //! Convergence criteria for linear solves
    double tolerance_{1.0e-8};

    //! Final relative residual norm of the linear system for the last iteration
    double finalResNorm_;

    //! Maximum number of iterations to perform when solving the system
    HYPRE_Int maxIterations_{200};

    //! Verbosity of HYPRE
    HYPRE_Int printLevel_{1};

    //! Logging verbosity of HYPRE
    HYPRE_Int logLevel_{1};

    //! Total number of iterations for convergence
    HYPRE_Int numIterations_;

    bool systemInitialized_{false};
    bool solverInitialized_{false};
    bool precondInitialized_{false};
    bool usePrecond_{false};

    bool isPeriodic_{true};
};

}  // nalu
}  // sierra


#endif /* HYPRESTRUCTSOLVER_H */
