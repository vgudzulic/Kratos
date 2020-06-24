//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//
//

// System includes
#include <functional>

// External includes

// Project includes
#include "containers/model.h"
#include "includes/cfd_variables.h"
#include "includes/model_part.h"
#include "includes/variables.h"

// Application includes
#include "custom_utilities/test_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "evm_k_epsilon_high_re_test_utilities.h"

namespace Kratos
{
namespace EvmKEpsilonHighReTestUtilities
{
ModelPart& RansEvmKEpsilonK2D3N_SetUp(Model& rModel, const std::string& rElementName)
{
    const auto add_variables_function = [](ModelPart& rModelPart) {
        rModelPart.AddNodalSolutionStepVariable(VELOCITY);
        rModelPart.AddNodalSolutionStepVariable(KINEMATIC_VISCOSITY);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_VISCOSITY);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_KINETIC_ENERGY);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_KINETIC_ENERGY_RATE);
        rModelPart.AddNodalSolutionStepVariable(RANS_AUXILIARY_VARIABLE_1);
    };

    using namespace RansModellingApplicationTestUtilities;

    ModelPart& r_model_part =
        CreateTestModelPart(rModel, rElementName, "LineCondition2D2N",
                            add_variables_function, TURBULENT_KINETIC_ENERGY, 1);

    // set nodal historical variables
    RandomFillNodalHistoricalVariable(r_model_part, VELOCITY, -10.0, 10.0);
    RandomFillNodalHistoricalVariable(r_model_part, KINEMATIC_VISCOSITY, 1e-5, 1e-3);
    RandomFillNodalHistoricalVariable(r_model_part, TURBULENT_VISCOSITY, 1e-3, 1e-1);
    RandomFillNodalHistoricalVariable(r_model_part, TURBULENT_KINETIC_ENERGY, 1.0, 100.0);
    RandomFillNodalHistoricalVariable(r_model_part, TURBULENT_KINETIC_ENERGY_RATE, 1.0, 50.0);
    RandomFillNodalHistoricalVariable(r_model_part, RANS_AUXILIARY_VARIABLE_1, 1.0, 10.0);

    // set process info variables
    ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
    r_process_info.SetValue(TURBULENT_KINETIC_ENERGY_SIGMA, 0.5);
    r_process_info.SetValue(TURBULENCE_RANS_C_MU, 2.1);

    return r_model_part;
}

ModelPart& RansEvmKEpsilonEpsilon2D3N_SetUp(Model& rModel, const std::string& rElementName)
{
    const auto add_variables_function = [](ModelPart& rModelPart) {
        rModelPart.AddNodalSolutionStepVariable(VELOCITY);
        rModelPart.AddNodalSolutionStepVariable(KINEMATIC_VISCOSITY);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_VISCOSITY);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_KINETIC_ENERGY);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_ENERGY_DISSIPATION_RATE);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_ENERGY_DISSIPATION_RATE_2);
        rModelPart.AddNodalSolutionStepVariable(RANS_AUXILIARY_VARIABLE_2);
    };

    using namespace RansModellingApplicationTestUtilities;

    ModelPart& r_model_part = CreateTestModelPart(
        rModel, rElementName, "LineCondition2D2N", add_variables_function,
        TURBULENT_ENERGY_DISSIPATION_RATE, 1);

    // set nodal historical variables
    RandomFillNodalHistoricalVariable(r_model_part, VELOCITY, -10.0, 10.0);
    RandomFillNodalHistoricalVariable(r_model_part, KINEMATIC_VISCOSITY, 1e-5, 1e-3);
    RandomFillNodalHistoricalVariable(r_model_part, TURBULENT_VISCOSITY, 1e-3, 1e-1);
    RandomFillNodalHistoricalVariable(r_model_part, TURBULENT_KINETIC_ENERGY, 1.0, 100.0);
    RandomFillNodalHistoricalVariable(
        r_model_part, TURBULENT_ENERGY_DISSIPATION_RATE, 1.0, 1000.0);
    RandomFillNodalHistoricalVariable(
        r_model_part, TURBULENT_ENERGY_DISSIPATION_RATE_2, 1.0, 1000.0);
    RandomFillNodalHistoricalVariable(r_model_part, RANS_AUXILIARY_VARIABLE_2, 1.0, 10.0);

    // set process info variables
    ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
    r_process_info.SetValue(TURBULENCE_RANS_C1, 3.1);
    r_process_info.SetValue(TURBULENCE_RANS_C2, 4.2);
    r_process_info.SetValue(TURBULENCE_RANS_C_MU, 1.1);
    r_process_info.SetValue(TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA, 1.1);

    return r_model_part;
}

ModelPart& RansEvmKEpsilonEpsilon2D2N_SetUp(Model& rModel, const std::string& rConditionName)
{
    const auto add_variables_function = [](ModelPart& rModelPart) {
        rModelPart.AddNodalSolutionStepVariable(VELOCITY);
        rModelPart.AddNodalSolutionStepVariable(KINEMATIC_VISCOSITY);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_VISCOSITY);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_KINETIC_ENERGY);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_ENERGY_DISSIPATION_RATE);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_ENERGY_DISSIPATION_RATE_2);
    };

    using namespace RansModellingApplicationTestUtilities;

    ModelPart& r_model_part =
        CreateTestModelPart(rModel, "Element2D3N", rConditionName, add_variables_function,
                            TURBULENT_ENERGY_DISSIPATION_RATE, 1);

    // set nodal historical variables
    RandomFillNodalHistoricalVariable(r_model_part, VELOCITY, -10.0, 10.0);
    RandomFillNodalHistoricalVariable(r_model_part, KINEMATIC_VISCOSITY, 1e-5, 1e-3);
    RandomFillNodalHistoricalVariable(r_model_part, TURBULENT_VISCOSITY, 1e-3, 1e-1);
    RandomFillNodalHistoricalVariable(r_model_part, TURBULENT_KINETIC_ENERGY, 10.0, 40.0);
    RandomFillNodalHistoricalVariable(
        r_model_part, TURBULENT_ENERGY_DISSIPATION_RATE, 1.0, 1000.0);
    RandomFillNodalHistoricalVariable(
        r_model_part, TURBULENT_ENERGY_DISSIPATION_RATE_2, 1.0, 1000.0);

    RandomFillContainerVariable<ModelPart::ConditionsContainerType>(
        r_model_part, RANS_Y_PLUS, 10.0, 100.0);

    // set process info variables
    ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
    r_process_info.SetValue(WALL_VON_KARMAN, 3.1);
    r_process_info.SetValue(WALL_SMOOTHNESS_BETA, 4.2);
    r_process_info.SetValue(RANS_Y_PLUS_LIMIT, 12.0);
    r_process_info.SetValue(TURBULENCE_RANS_C_MU, 0.09);
    r_process_info.SetValue(TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA, 1.1);

    return r_model_part;
}
} // namespace EvmKEpsilonHighReTestUtilities
} // namespace Kratos
