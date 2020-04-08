//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

// System includes

// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "custom_python/add_custom_strategies_to_python.h"
#include "includes/define_python.h"
#include "linear_solvers/linear_solver.h"
#include "spaces/ublas_space.h"

// strategies
#include "custom_strategies/coupled_strategy.h"
#include "custom_strategies/coupled_strategy_item.h"

// schemes
#include "custom_strategies/generic_residual_based_bossak_velocity_scalar_scheme.h"
#include "custom_strategies/generic_residualbased_simple_steady_scalar_scheme.h"
#include "custom_strategies/residualbased_simple_steady_velocity_scheme.h"
#include "custom_strategies/algebraic_flux_corrected_scalar_steady_scheme.h"

// convergence criterians
#include "custom_strategies/generic_convergence_criteria.h"

namespace Kratos
{
namespace Python
{
void AddCustomStrategiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    using LocalSpaceType = UblasSpace<double, Matrix, Vector>;
    using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
    using BaseSchemeType = Scheme<SparseSpaceType, LocalSpaceType>;
    using LinearSolverType = LinearSolver<SparseSpaceType, LocalSpaceType>;
    using BaseSolvingStrategyType =
        SolvingStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>;

    // Add coupled strategy item
    using CoupledStrategyItemType =
        CoupledStrategyItem<SparseSpaceType, LocalSpaceType, LinearSolverType>;
    py::class_<CoupledStrategyItemType, typename CoupledStrategyItemType::Pointer>(
        m, "CoupledStrategyItem")
        .def(py::init<BaseSolvingStrategyType::Pointer, std::string, int>())
        .def(py::init<BaseSolvingStrategyType::Pointer, std::string, int, std::vector<int>>())
        .def("AddAuxiliaryProcess", &CoupledStrategyItemType::AddAuxiliaryProcess)
        .def("GetName", &CoupledStrategyItemType::GetName)
        .def("GetStrategy", &CoupledStrategyItemType::GetStrategy)
        .def("GetAuxiliaryProcessList", &CoupledStrategyItemType::GetStrategy)
        .def("GetStrategyInfo", &CoupledStrategyItemType::GetStrategyInfo)
        .def("GetStrategySolvabilityPattern", &CoupledStrategyItemType::GetStrategySolvabilityPattern)
        .def("SetStrategySolvabilityPattern", &CoupledStrategyItemType::SetStrategySolvabilityPattern);

    // Add strtegies
    using CoupledStrategyType =
        CoupledStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>;
    py::class_<CoupledStrategyType, typename CoupledStrategyType::Pointer, BaseSolvingStrategyType>(
        m, "CoupledStrategy")
        .def(py::init<ModelPart&, bool, bool, bool, int>())
        .def("AddStrategyItem", &CoupledStrategyType::AddStrategyItem)
        .def("AddConvergenceCheckVariable", &CoupledStrategyType::AddConvergenceCheckVariable);

    // Convergence criteria
    py::class_<GenericConvergenceCriteria<SparseSpaceType, LocalSpaceType>,
               typename GenericConvergenceCriteria<SparseSpaceType, LocalSpaceType>::Pointer,
               ConvergenceCriteria<SparseSpaceType, LocalSpaceType>>(
        m, "GenericScalarConvergenceCriteria")
        .def(py::init<double, double>());

    py::class_<GenericResidualBasedBossakVelocityScalarScheme<SparseSpaceType, LocalSpaceType>,
               typename GenericResidualBasedBossakVelocityScalarScheme<SparseSpaceType, LocalSpaceType>::Pointer, BaseSchemeType>(
        m, "GenericResidualBasedBossakVelocityDynamicScalarScheme")
        .def(py::init<const double, const double, const Variable<double>&,
                      const Variable<double>&, const Variable<double>&>());

    py::class_<GenericResidualBasedSimpleSteadyScalarScheme<SparseSpaceType, LocalSpaceType>,
               typename GenericResidualBasedSimpleSteadyScalarScheme<SparseSpaceType, LocalSpaceType>::Pointer, BaseSchemeType>(
        m, "GenericResidualBasedSimpleSteadyScalarScheme")
        .def(py::init<const double>());

    py::class_<ResidualBasedSimpleSteadyVelocityScheme<SparseSpaceType, LocalSpaceType>,
               typename ResidualBasedSimpleSteadyVelocityScheme<SparseSpaceType, LocalSpaceType>::Pointer, BaseSchemeType>(
        m, "ResidualBasedSimpleSteadyVelocityScheme")
        .def(py::init<const double, const unsigned int>());

    py::class_<AlgebraicFluxCorrectedScalarSteadyScheme<SparseSpaceType, LocalSpaceType>,
               typename AlgebraicFluxCorrectedScalarSteadyScheme<SparseSpaceType, LocalSpaceType>::Pointer, BaseSchemeType>(
        m, "AlgebraicFluxCorrectedScalarSteadyScheme")
        .def(py::init<const double>());
}

} // namespace Python.
} // Namespace Kratos
