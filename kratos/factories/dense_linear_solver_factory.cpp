//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//			 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// Project includes
#include "includes/define.h"
#include "linear_solvers/linear_solver.h"

// Linear solvers
#include "factories/dense_linear_solver_factory.h"
#include "linear_solvers/dense_lu_direct_solver.h"

namespace Kratos {

void RegisterDenseLinearSolvers()
{
    typedef TUblasDenseSpace<double> LocalSpaceType;

    // Dense LU solver
    static auto DenseLUDirectSolverFactory = DenseLUDirectSolver<LocalSpaceType, LocalSpaceType>::Factory();
    KRATOS_REGISTER_DENSE_LINEAR_SOLVER("dense_lu_solver", DenseLUDirectSolverFactory);

}

template class KratosComponents<DenseLinearSolverFactoryType>;
template class KratosComponents<ComplexDenseLinearSolverFactoryType>;
}
