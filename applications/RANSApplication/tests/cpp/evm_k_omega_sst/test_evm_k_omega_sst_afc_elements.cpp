//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:
//

// System includes
#include <string>
#include <vector>

// External includes

// Project includes
#include "containers/model.h"
#include "includes/cfd_variables.h"
#include "includes/kratos_components.h"
#include "includes/variables.h"
#include "testing/testing.h"

// Application includes
#include "custom_utilities/test_utilities.h"
#include "rans_application_variables.h"

namespace Kratos
{
namespace Testing
{
namespace
{
ModelPart& RansEvmKOmegaSSTKAFC2D3N_SetUp(Model& rModel)
{
    const auto add_variables_function = [](ModelPart& rModelPart) {
        rModelPart.AddNodalSolutionStepVariable(VELOCITY);
        rModelPart.AddNodalSolutionStepVariable(KINEMATIC_VISCOSITY);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_VISCOSITY);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_KINETIC_ENERGY);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_KINETIC_ENERGY_RATE);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE);
        rModelPart.AddNodalSolutionStepVariable(RANS_AUXILIARY_VARIABLE_1);
        rModelPart.AddNodalSolutionStepVariable(DISTANCE);
    };

    using namespace RansModellingApplicationTestUtilities;

    ModelPart& r_model_part = CreateTestModelPart(
        rModel, "RansEvmKOmegaSSTKAFC2D3N", "LineCondition2D2N",
        add_variables_function, TURBULENT_KINETIC_ENERGY, 1);

    // set nodal historical variables
    RandomFillNodalHistoricalVariable(r_model_part, VELOCITY, -10.0, 10.0);
    RandomFillNodalHistoricalVariable(r_model_part, KINEMATIC_VISCOSITY, 1e-5, 1e-3);
    RandomFillNodalHistoricalVariable(r_model_part, TURBULENT_VISCOSITY, 1e-3, 1e-1);
    RandomFillNodalHistoricalVariable(r_model_part, TURBULENT_KINETIC_ENERGY, 1.0, 100.0);
    RandomFillNodalHistoricalVariable(r_model_part, TURBULENT_KINETIC_ENERGY_RATE, 1.0, 50.0);
    RandomFillNodalHistoricalVariable(
        r_model_part, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, 1.0, 1000.0);
    RandomFillNodalHistoricalVariable(r_model_part, RANS_AUXILIARY_VARIABLE_1, 1.0, 10.0);
    RandomFillNodalHistoricalVariable(r_model_part, DISTANCE, 1.0, 6.0);

    // set process info variables
    ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
    r_process_info.SetValue(TURBULENT_KINETIC_ENERGY_SIGMA_1, 0.5);
    r_process_info.SetValue(TURBULENT_KINETIC_ENERGY_SIGMA_2, 0.3);
    r_process_info.SetValue(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA_2, 2.0);
    r_process_info.SetValue(TURBULENCE_RANS_C_MU, 2.1);

    return r_model_part;
}

ModelPart& RansEvmKOmegaSSTOmegaAFC2D3N_SetUp(Model& rModel)
{
    const auto add_variables_function = [](ModelPart& rModelPart) {
        rModelPart.AddNodalSolutionStepVariable(VELOCITY);
        rModelPart.AddNodalSolutionStepVariable(KINEMATIC_VISCOSITY);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_VISCOSITY);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_KINETIC_ENERGY);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_2);
        rModelPart.AddNodalSolutionStepVariable(RANS_AUXILIARY_VARIABLE_2);
        rModelPart.AddNodalSolutionStepVariable(DISTANCE);
    };

    using namespace RansModellingApplicationTestUtilities;

    ModelPart& r_model_part = CreateTestModelPart(
        rModel, "RansEvmKOmegaSSTOmegaAFC2D3N", "LineCondition2D2N",
        add_variables_function, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, 1);

    // set nodal historical variables
    RandomFillNodalHistoricalVariable(r_model_part, VELOCITY, -10.0, 10.0);
    RandomFillNodalHistoricalVariable(r_model_part, KINEMATIC_VISCOSITY, 1e-5, 1e-3);
    RandomFillNodalHistoricalVariable(r_model_part, TURBULENT_VISCOSITY, 1e-3, 1e-1);
    RandomFillNodalHistoricalVariable(r_model_part, TURBULENT_KINETIC_ENERGY, 1.0, 100.0);
    RandomFillNodalHistoricalVariable(
        r_model_part, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, 1.0, 1000.0);
    RandomFillNodalHistoricalVariable(
        r_model_part, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_2, 1.0, 1000.0);
    RandomFillNodalHistoricalVariable(r_model_part, RANS_AUXILIARY_VARIABLE_2, 1.0, 10.0);
    RandomFillNodalHistoricalVariable(r_model_part, DISTANCE, 1.0, 6.0);

    // set process info variables
    ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
    r_process_info.SetValue(TURBULENCE_RANS_BETA_1, 3.1);
    r_process_info.SetValue(TURBULENCE_RANS_BETA_2, 4.2);
    r_process_info.SetValue(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA_1, 1.1);
    r_process_info.SetValue(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA_2, 2.1);
    r_process_info.SetValue(TURBULENCE_RANS_C_MU, 0.4);
    r_process_info.SetValue(WALL_VON_KARMAN, 5.2);

    return r_model_part;
}

} // namespace

KRATOS_TEST_CASE_IN_SUITE(RansEvmKOmegaSSTKAFC2D3N_EquationIdVector, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansEvmKOmegaSSTKAFC2D3N_SetUp(model);

    // Test:
    RansModellingApplicationTestUtilities::TestEquationIdVector(r_model_part);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKOmegaSSTKAFC2D3N_GetDofList, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansEvmKOmegaSSTKAFC2D3N_SetUp(model);

    // Test:
    RansModellingApplicationTestUtilities::TestGetDofList(r_model_part, TURBULENT_KINETIC_ENERGY);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKOmegaSSTKAFC2D3N_CalculateLocalSystem, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansEvmKOmegaSSTKAFC2D3N_SetUp(model);

    // Test:
    Matrix LHS, ref_LHS(3, 3);
    Vector RHS, ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateLocalSystem(LHS, RHS, r_model_part.GetProcessInfo());

    // setting reference values
    ref_RHS[0] = 2.67829732581686741355e+00;
    ref_RHS[1] = 2.21576728463429173388e+00;
    ref_RHS[2] = 1.85563938248039850265e+00;
    ref_LHS = ZeroMatrix(3, 3);

    std::cout << std::scientific << std::setprecision(20) << std::endl
              << RHS[0] << std::endl
              << RHS[1] << std::endl
              << RHS[2] << std::endl;

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKOmegaSSTKAFC2D3N_CalculateRightHandSide, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansEvmKOmegaSSTKAFC2D3N_SetUp(model);

    // Test:
    Vector RHS, ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateRightHandSide(RHS, r_model_part.GetProcessInfo());

    // setting reference values
    ref_RHS[0] = 2.67829732581686741355e+00;
    ref_RHS[1] = 2.21576728463429173388e+00;
    ref_RHS[2] = 1.85563938248039850265e+00;

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKOmegaSSTKAFC2D3N_CalculateLocalVelocityContribution, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansEvmKOmegaSSTKAFC2D3N_SetUp(model);

    // Test:
    Matrix LHS, ref_LHS(3, 3);
    Vector RHS, ref_RHS;
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateLocalVelocityContribution(LHS, RHS, r_model_part.GetProcessInfo());

    // setting reference values
    ref_LHS(0, 0) = 1.22721919178327041777e+02;
    ref_LHS(0, 1) = 7.87547318032734722237e+01;
    ref_LHS(0, 2) = 8.08415536207580629480e+01;
    ref_LHS(1, 0) = 7.99174392233141759334e+01;
    ref_LHS(1, 1) = 1.89298008826390287140e+02;
    ref_LHS(1, 2) = 9.45742683315857988191e+01;
    ref_LHS(2, 0) = 8.11121632251401223357e+01;
    ref_LHS(2, 1) = 9.44390118548673456189e+01;
    ref_LHS(2, 2) = 1.98028509313836764250e+02;

    std::cout << std::scientific << std::setprecision(20) << std::endl
              << LHS(0, 0) << std::endl
              << LHS(0, 1) << std::endl
              << LHS(0, 2) << std::endl
              << LHS(1, 0) << std::endl
              << LHS(1, 1) << std::endl
              << LHS(1, 2) << std::endl
              << LHS(2, 0) << std::endl
              << LHS(2, 1) << std::endl
              << LHS(2, 2) << std::endl;

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKOmegaSSTKAFC2D3N_CalculateMassMatrix, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansEvmKOmegaSSTKAFC2D3N_SetUp(model);

    // Test:
    Matrix M, ref_M;
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateMassMatrix(M, r_model_part.GetProcessInfo());

    // setting reference values
    ref_M = ZeroMatrix(3, 3);
    ref_M(0, 0) = 1.66666666666666657415e-01;
    ref_M(1, 1) = 1.66666666666666657415e-01;
    ref_M(2, 2) = 1.66666666666666657415e-01;

    std::cout << std::scientific << std::setprecision(20) << std::endl
              << M(0, 0) << std::endl
              << M(1, 1) << std::endl
              << M(2, 2) << std::endl;

    KRATOS_CHECK_MATRIX_NEAR(M, ref_M, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKOmegaSSTKAFC2D3N_CalculateDampingMatrix, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansEvmKOmegaSSTKAFC2D3N_SetUp(model);

    // Test:
    Matrix D, ref_D(3, 3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateDampingMatrix(D, r_model_part.GetProcessInfo());

    // setting reference values
    ref_D(0, 0) = 1.22721919178327041777e+02;
    ref_D(0, 1) = 7.87547318032734722237e+01;
    ref_D(0, 2) = 8.08415536207580629480e+01;
    ref_D(1, 0) = 7.99174392233141759334e+01;
    ref_D(1, 1) = 1.89298008826390287140e+02;
    ref_D(1, 2) = 9.45742683315857988191e+01;
    ref_D(2, 0) = 8.11121632251401223357e+01;
    ref_D(2, 1) = 9.44390118548673456189e+01;
    ref_D(2, 2) = 1.98028509313836764250e+02;

    KRATOS_CHECK_MATRIX_NEAR(D, ref_D, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKOmegaSSTOmegaAFC2D3N_EquationIdVector, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansEvmKOmegaSSTOmegaAFC2D3N_SetUp(model);

    // Test:
    RansModellingApplicationTestUtilities::TestEquationIdVector(r_model_part);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKOmegaSSTOmegaAFC2D3N_GetDofList, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansEvmKOmegaSSTOmegaAFC2D3N_SetUp(model);

    // Test:
    RansModellingApplicationTestUtilities::TestGetDofList(
        r_model_part, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKOmegaSSTOmegaAFC2D3N_CalculateLocalSystem, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansEvmKOmegaSSTOmegaAFC2D3N_SetUp(model);

    // Test:
    Matrix LHS, ref_LHS(3, 3);
    Vector RHS, ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateLocalSystem(LHS, RHS, r_model_part.GetProcessInfo());

    // setting reference values
    ref_RHS[0] = -3.16970054379328485084e+03;
    ref_RHS[1] = -3.16970054861076550878e+03;
    ref_RHS[2] = -3.16970020196045061311e+03;

    ref_LHS = ZeroMatrix(3, 3);

    std::cout << std::scientific << std::setprecision(20) << std::endl
              << RHS[0] << std::endl
              << RHS[1] << std::endl
              << RHS[2] << std::endl;

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKOmegaSSTOmegaAFC2D3N_CalculateRightHandSide, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansEvmKOmegaSSTOmegaAFC2D3N_SetUp(model);

    // Test:
    Vector RHS, ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateRightHandSide(RHS, r_model_part.GetProcessInfo());

    // setting reference values
    ref_RHS[0] = -3.16970054379328485084e+03;
    ref_RHS[1] = -3.16970054861076550878e+03;
    ref_RHS[2] = -3.16970020196045061311e+03;

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKOmegaSSTOmegaAFC2D3N_CalculateLocalVelocityContribution,
                          KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansEvmKOmegaSSTOmegaAFC2D3N_SetUp(model);

    // Test:
    Matrix LHS, ref_LHS(3, 3);
    Vector RHS, ref_RHS;
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateLocalVelocityContribution(LHS, RHS, r_model_part.GetProcessInfo());

    // setting reference values
    ref_LHS(0, 0) = 1.20138773411655847667e+02;
    ref_LHS(0, 1) = 7.13517563520743323124e+01;
    ref_LHS(0, 2) = 5.04570356145520690916e+01;
    ref_LHS(1, 0) = 7.25144637721150502330e+01;
    ref_LHS(1, 1) = 1.78527771954713898594e+02;
    ref_LHS(1, 2) = 6.24915534404995440809e+01;
    ref_LHS(2, 0) = 5.07276452189341497956e+01;
    ref_LHS(2, 1) = 6.23562969637810837753e+01;
    ref_LHS(2, 2) = 7.20466663477299107399e+01;

    std::cout << std::scientific << std::setprecision(20) << std::endl
              << LHS(0, 0) << std::endl
              << LHS(0, 1) << std::endl
              << LHS(0, 2) << std::endl
              << LHS(1, 0) << std::endl
              << LHS(1, 1) << std::endl
              << LHS(1, 2) << std::endl
              << LHS(2, 0) << std::endl
              << LHS(2, 1) << std::endl
              << LHS(2, 2) << std::endl;

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKOmegaSSTOmegaAFC2D3N_CalculateMassMatrix, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansEvmKOmegaSSTOmegaAFC2D3N_SetUp(model);

    // Test:
    Matrix M, ref_M;
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateMassMatrix(M, r_model_part.GetProcessInfo());

    // setting reference values
    ref_M = ZeroMatrix(3, 3);
    ref_M(0, 0) = 1.66666666666666657415e-01;
    ref_M(1, 1) = 1.66666666666666657415e-01;
    ref_M(2, 2) = 1.66666666666666657415e-01;

    std::cout << std::scientific << std::setprecision(20) << std::endl
              << M(0, 0) << std::endl
              << M(1, 1) << std::endl
              << M(2, 2) << std::endl;

    KRATOS_CHECK_MATRIX_NEAR(M, ref_M, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKOmegaSSTOmegaAFC2D3N_CalculateDampingMatrix, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansEvmKOmegaSSTOmegaAFC2D3N_SetUp(model);

    // Test:
    Matrix D, ref_D(3, 3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateDampingMatrix(D, r_model_part.GetProcessInfo());

    // setting reference values
    ref_D(0, 0) = 1.20138773411655847667e+02;
    ref_D(0, 1) = 7.13517563520743323124e+01;
    ref_D(0, 2) = 5.04570356145520690916e+01;
    ref_D(1, 0) = 7.25144637721150502330e+01;
    ref_D(1, 1) = 1.78527771954713898594e+02;
    ref_D(1, 2) = 6.24915534404995440809e+01;
    ref_D(2, 0) = 5.07276452189341497956e+01;
    ref_D(2, 1) = 6.23562969637810837753e+01;
    ref_D(2, 2) = 7.20466663477299107399e+01;

    KRATOS_CHECK_MATRIX_NEAR(D, ref_D, 1e-12);
}

} // namespace Testing
} // namespace Kratos.
