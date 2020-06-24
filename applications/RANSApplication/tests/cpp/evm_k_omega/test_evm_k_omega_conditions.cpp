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

// External includes

// Project includes
#include "containers/model.h"
#include "testing/testing.h"

// Application includes
#include "custom_utilities/test_utilities.h"
#include "evm_k_omega_test_utilities.h"
#include "rans_application_variables.h"

namespace Kratos
{
namespace Testing
{
namespace
{
ModelPart& RansEvmKOmegaOmegaKBasedWallCondition2D2N_SetUp(Model& rModel)
{
    ModelPart& r_model_part = EvmKOmegaTestUtilities::RansEvmKOmegaOmega2D2N_SetUp(
        rModel, "RansEvmKOmegaOmegaKBasedWallCondition2D2N");

    return r_model_part;
}

ModelPart& RansEvmKOmegaOmegaUBasedWallCondition2D2N_SetUp(Model& rModel)
{
    ModelPart& r_model_part = EvmKOmegaTestUtilities::RansEvmKOmegaOmega2D2N_SetUp(
        rModel, "RansEvmKOmegaOmegaUBasedWallCondition2D2N");

    return r_model_part;
}

} // namespace

KRATOS_TEST_CASE_IN_SUITE(RansEvmKOmegaOmegaKBasedWallCondition2D2N_EquationIdVector, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansEvmKOmegaOmegaKBasedWallCondition2D2N_SetUp(model);

    // Test:
    RansModellingApplicationTestUtilities::TestEquationIdVector<ModelPart::ConditionsContainerType>(
        r_model_part);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKOmegaOmegaKBasedWallCondition2D2N_GetDofList, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansEvmKOmegaOmegaKBasedWallCondition2D2N_SetUp(model);

    // Test:
    RansModellingApplicationTestUtilities::TestGetDofList<ModelPart::ConditionsContainerType>(
        r_model_part, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKOmegaOmegaKBasedWallCondition2D2N_CalculateLocalSystem,
                          KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansEvmKOmegaOmegaKBasedWallCondition2D2N_SetUp(model);

    // Test:
    Matrix LHS, ref_LHS;
    Vector RHS, ref_RHS;
    auto& r_condition = r_model_part.Conditions().front();

    // checking for no-wall function
    r_condition.SetValue(RANS_IS_WALL, 0);
    r_condition.CalculateLocalSystem(LHS, RHS, r_model_part.GetProcessInfo());
    // setting reference values
    ref_RHS = ZeroVector(2);
    ref_LHS = ZeroMatrix(2, 2);

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);

    // checking for wall function
    r_condition.SetValue(RANS_IS_WALL, 1);
    r_condition.CalculateLocalSystem(LHS, RHS, r_model_part.GetProcessInfo());
    // setting reference values
    ref_RHS[0] = 1.9460396765785567e+03;
    ref_RHS[1] = 1.9460396765785567e+03;
    ref_LHS = ZeroMatrix(2, 2);

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKOmegaOmegaKBasedWallCondition2D2N_CalculateRightHandSide,
                          KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansEvmKOmegaOmegaKBasedWallCondition2D2N_SetUp(model);

    // Test:
    Vector RHS, ref_RHS;
    auto& r_condition = r_model_part.Conditions().front();

    // checking for no-wall function
    r_condition.SetValue(RANS_IS_WALL, 0);
    r_condition.CalculateRightHandSide(RHS, r_model_part.GetProcessInfo());
    // setting reference values
    ref_RHS = ZeroVector(2);

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);

    // checking for wall function
    r_condition.SetValue(RANS_IS_WALL, 1);
    r_condition.CalculateRightHandSide(RHS, r_model_part.GetProcessInfo());
    // setting reference values
    ref_RHS[0] = 1.9460396765785567e+03;
    ref_RHS[1] = 1.9460396765785567e+03;

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKOmegaOmegaUBasedWallCondition2D2N_EquationIdVector, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansEvmKOmegaOmegaUBasedWallCondition2D2N_SetUp(model);

    // Test:
    RansModellingApplicationTestUtilities::TestEquationIdVector<ModelPart::ConditionsContainerType>(
        r_model_part);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKOmegaOmegaUBasedWallCondition2D2N_GetDofList, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansEvmKOmegaOmegaUBasedWallCondition2D2N_SetUp(model);

    // Test:
    RansModellingApplicationTestUtilities::TestGetDofList<ModelPart::ConditionsContainerType>(
        r_model_part, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKOmegaOmegaUBasedWallCondition2D2N_CalculateLocalSystem,
                          KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansEvmKOmegaOmegaUBasedWallCondition2D2N_SetUp(model);

    // Test:
    Matrix LHS, ref_LHS;
    Vector RHS, ref_RHS;
    auto& r_condition = r_model_part.Conditions().front();

    // checking for no-wall function
    r_condition.SetValue(RANS_IS_WALL, 0);
    r_condition.CalculateLocalSystem(LHS, RHS, r_model_part.GetProcessInfo());
    // setting reference values
    ref_RHS = ZeroVector(2);
    ref_LHS = ZeroMatrix(2, 2);

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);

    // checking for wall function
    r_condition.SetValue(RANS_IS_WALL, 1);
    r_condition.CalculateLocalSystem(LHS, RHS, r_model_part.GetProcessInfo());
    // setting reference values
    ref_RHS[0] = 8.3262170750205513e+01;
    ref_RHS[1] = 8.3262170750205513e+01;
    ref_LHS = ZeroMatrix(2, 2);

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKOmegaOmegaUBasedWallCondition2D2N_CalculateRightHandSide,
                          KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansEvmKOmegaOmegaUBasedWallCondition2D2N_SetUp(model);

    // Test:
    Vector RHS, ref_RHS;
    auto& r_condition = r_model_part.Conditions().front();

    // checking for no-wall function
    r_condition.SetValue(RANS_IS_WALL, 0);
    r_condition.CalculateRightHandSide(RHS, r_model_part.GetProcessInfo());
    // setting reference values
    ref_RHS = ZeroVector(2);

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);

    // checking for wall function
    r_condition.SetValue(RANS_IS_WALL, 1);
    r_condition.CalculateRightHandSide(RHS, r_model_part.GetProcessInfo());
    // setting reference values
    ref_RHS[0] = 8.3262170750205513e+01;
    ref_RHS[1] = 8.3262170750205513e+01;

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
}

} // namespace Testing
} // namespace Kratos.
