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
#include <vector>

// External includes

// Project includes
#include "containers/model.h"
#include "testing/testing.h"

// Application includes
#include "custom_elements/evm_k_epsilon/evm_k_epsilon_adjoint_utilities.h"
#include "custom_processes/auxiliary_processes/rans_logarithmic_y_plus_calculation_process.h"
#include "custom_processes/auxiliary_processes/rans_logarithmic_y_plus_velocity_sensitivities_process.h"
#include "custom_processes/auxiliary_processes/rans_nut_k_epsilon_high_re_calculation_process.h"
#include "custom_processes/auxiliary_processes/rans_nut_k_epsilon_high_re_sensitivities_process.h"
#include "custom_utilities/rans_calculation_utilities.h"
#include "custom_utilities/test_utilities.h"
#include "test_k_epsilon_utilities.h"

namespace Kratos
{
namespace Testing
{
KRATOS_TEST_CASE_IN_SUITE(RansEvmEpsilonAdjoint2D3N_EquationIdVector, RANSEvModelsKEpsilonElementMethods)
{
    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonElementTestModelPart(
        r_adjoint_model_part, "RansEvmEpsilonAdjoint2D3N");

    for (auto& element : r_adjoint_model_part.Elements())
    {
        std::vector<std::size_t> equation_ids{};
        element.EquationIdVector(equation_ids, r_adjoint_model_part.GetProcessInfo());
        KRATOS_CHECK_EQUAL(equation_ids.size(), element.GetGeometry().PointsNumber());

        for (std::size_t i = 0; i < equation_ids.size(); ++i)
        {
            KRATOS_ERROR_IF(
                equation_ids[i] !=
                element.GetGeometry()[i].GetDof(RANS_SCALAR_2_ADJOINT_1).EquationId())
                << "Equation id mismatch.";
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmEpsilonAdjoint2D3N_GetDofList, RANSEvModelsKEpsilonElementMethods)
{
    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonElementTestModelPart(
        r_adjoint_model_part, "RansEvmEpsilonAdjoint2D3N");

    for (auto& element : r_adjoint_model_part.Elements())
    {
        auto dofs = Element::DofsVectorType{};
        element.GetDofList(dofs, r_adjoint_model_part.GetProcessInfo());
        KRATOS_CHECK_EQUAL(dofs.size(), element.GetGeometry().PointsNumber());
        for (std::size_t i = 0; i < dofs.size(); ++i)
        {
            KRATOS_ERROR_IF(
                dofs[i] !=
                element.GetGeometry()[i].pGetDof(RANS_SCALAR_2_ADJOINT_1))
                << "Dofs mismatch.";
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmEpsilonAdjoint2D3N_GetValuesVector,
                          RANSEvModelsKEpsilonElementMethods)
{
    Model adjoint_model;
    ModelPart& r_adjoint_model_part =
        adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonElementTestModelPart(
        r_adjoint_model_part, "RansEvmEpsilonAdjoint2D3N");

    for (IndexType i_element = 0;
         i_element < r_adjoint_model_part.NumberOfElements(); ++i_element)
    {
        Element& r_element = *(r_adjoint_model_part.ElementsBegin() + i_element);
        GeometryType& r_geometry = r_element.GetGeometry();
        const IndexType number_of_nodes = r_geometry.PointsNumber();

        Vector element_values;
        r_element.GetValuesVector(element_values);

        Vector values(3);
        IndexType local_index = 0;
        for (IndexType i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            const NodeType& r_node = r_geometry[i_node];
            values[local_index++] = r_node.FastGetSolutionStepValue(RANS_SCALAR_2_ADJOINT_1);
        }

        RansModellingApplicationTestUtilities::CheckNear(element_values, values);
    }
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmEpsilonAdjoint2D3N_GetFirstDerivativesVector,
                          RANSEvModelsKEpsilonElementMethods)
{
    Model adjoint_model;
    ModelPart& r_adjoint_model_part =
        adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonElementTestModelPart(
        r_adjoint_model_part, "RansEvmEpsilonAdjoint2D3N");

    for (IndexType i_element = 0;
         i_element < r_adjoint_model_part.NumberOfElements(); ++i_element)
    {
        Element& r_element = *(r_adjoint_model_part.ElementsBegin() + i_element);

        Vector element_values;
        r_element.GetFirstDerivativesVector(element_values);

        Vector values = ZeroVector(3);

        RansModellingApplicationTestUtilities::CheckNear(element_values, values);
    }
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmEpsilonAdjoint2D3N_GetSecondDerivativesVector,
                          RANSEvModelsKEpsilonElementMethods)
{
    Model adjoint_model;
    ModelPart& r_adjoint_model_part =
        adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonElementTestModelPart(
        r_adjoint_model_part, "RansEvmEpsilonAdjoint2D3N");

    for (IndexType i_element = 0;
         i_element < r_adjoint_model_part.NumberOfElements(); ++i_element)
    {
        Element& r_element = *(r_adjoint_model_part.ElementsBegin() + i_element);
        GeometryType& r_geometry = r_element.GetGeometry();
        const IndexType number_of_nodes = r_geometry.PointsNumber();

        Vector element_values;
        r_element.GetSecondDerivativesVector(element_values);

        Vector values(3);
        IndexType local_index = 0;
        for (IndexType i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            const NodeType& r_node = r_geometry[i_node];
            values[local_index++] =
                r_node.FastGetSolutionStepValue(RANS_SCALAR_2_ADJOINT_3);
        }

        RansModellingApplicationTestUtilities::CheckNear(element_values, values);
    }
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmEpsilonAdjoint2D3N_CalculateFirstDerivativesLHS,
                          RANSModellingApplicationElementInterfaces)
{
    Model primal_model;
    ModelPart& r_primal_model_part = primal_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonElementTestModelPart(
        r_primal_model_part, "RansEvmEpsilon2D3N");

    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonElementTestModelPart(
        r_adjoint_model_part, "RansEvmEpsilonAdjoint2D3N");

    Parameters empty_y_plus_parameters = Parameters(R"({
        "model_part_name" : "test"
    })");
    RansLogarithmicYPlusVelocitySensitivitiesProcess y_plus_sensitivities_process(
        adjoint_model, empty_y_plus_parameters);
    RansLogarithmicYPlusCalculationProcess adjoint_y_plus_process(
        adjoint_model, empty_y_plus_parameters);
    RansLogarithmicYPlusCalculationProcess primal_y_plus_process(
        primal_model, empty_y_plus_parameters);

    Parameters empty_nut_parameters = Parameters(R"({
        "model_part_name" : "test"
    })");
    RansNutKEpsilonHighReSensitivitiesProcess nut_sensitivities_process(
        adjoint_model, empty_nut_parameters);
    RansNutKEpsilonHighReCalculationProcess adjoint_nut_process(
        adjoint_model, empty_nut_parameters);
    RansNutKEpsilonHighReCalculationProcess primal_nut_process(primal_model, empty_nut_parameters);

    auto perturb_variable = [](NodeType& rNode) -> double& {
        return rNode.FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE);
    };

    auto calculate_sensitivity_matrix = [](Matrix& rOutput, Element& rElement,
                                           ProcessInfo& rProcessInfo) {
        rElement.CalculateFirstDerivativesLHS(rOutput, rProcessInfo);
    };

    RansModellingApplicationTestUtilities::RunResidualScalarSensitivityTest(
        r_primal_model_part, r_adjoint_model_part, primal_y_plus_process, primal_nut_process,
        adjoint_y_plus_process, adjoint_nut_process, y_plus_sensitivities_process,
        nut_sensitivities_process, RansEvmKEpsilonModel::UpdateVariablesInModelPart,
        calculate_sensitivity_matrix, perturb_variable, 1e-7, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmEpsilonAdjoint2D3N_Calculate_RANS_TURBULENT_KINETIC_ENERGY_PARTIAL_DERIVATIVE,
                          RANSModellingApplicationElementInterfaces)
{
    Model primal_model;
    ModelPart& r_primal_model_part = primal_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonElementTestModelPart(
        r_primal_model_part, "RansEvmEpsilon2D3N");

    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonElementTestModelPart(
        r_adjoint_model_part, "RansEvmEpsilonAdjoint2D3N");

    Parameters empty_y_plus_parameters = Parameters(R"({
        "model_part_name" : "test"
    })");
    RansLogarithmicYPlusVelocitySensitivitiesProcess y_plus_sensitivities_process(
        adjoint_model, empty_y_plus_parameters);
    RansLogarithmicYPlusCalculationProcess adjoint_y_plus_process(
        adjoint_model, empty_y_plus_parameters);
    RansLogarithmicYPlusCalculationProcess primal_y_plus_process(
        primal_model, empty_y_plus_parameters);

    Parameters empty_nut_parameters = Parameters(R"({
        "model_part_name" : "test"
    })");
    RansNutKEpsilonHighReSensitivitiesProcess nut_sensitivities_process(
        adjoint_model, empty_nut_parameters);
    RansNutKEpsilonHighReCalculationProcess adjoint_nut_process(
        adjoint_model, empty_nut_parameters);
    RansNutKEpsilonHighReCalculationProcess primal_nut_process(primal_model, empty_nut_parameters);

    auto perturb_variable = [](NodeType& rNode) -> double& {
        return rNode.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY);
    };

    auto calculate_sensitivity_matrix = [](Matrix& rOutput, Element& rElement,
                                           ProcessInfo& rProcessInfo) {
        rElement.Calculate(RANS_TURBULENT_KINETIC_ENERGY_PARTIAL_DERIVATIVE,
                           rOutput, rProcessInfo);
    };

    RansModellingApplicationTestUtilities::RunResidualScalarSensitivityTest(
        r_primal_model_part, r_adjoint_model_part, primal_y_plus_process, primal_nut_process,
        adjoint_y_plus_process, adjoint_nut_process, y_plus_sensitivities_process,
        nut_sensitivities_process, RansEvmKEpsilonModel::UpdateVariablesInModelPart,
        calculate_sensitivity_matrix, perturb_variable, 1e-8, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmEpsilonAdjoint2D3N_CalculateSecondDerivativesLHS,
                          RANSModellingApplicationElementInterfaces)
{
    Model primal_model;
    ModelPart& r_primal_model_part = primal_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonElementTestModelPart(
        r_primal_model_part, "RansEvmEpsilon2D3N");

    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonElementTestModelPart(
        r_adjoint_model_part, "RansEvmEpsilonAdjoint2D3N");

    Parameters empty_y_plus_parameters = Parameters(R"({
        "model_part_name" : "test"
    })");
    RansLogarithmicYPlusVelocitySensitivitiesProcess y_plus_sensitivities_process(
        adjoint_model, empty_y_plus_parameters);
    RansLogarithmicYPlusCalculationProcess adjoint_y_plus_process(
        adjoint_model, empty_y_plus_parameters);
    RansLogarithmicYPlusCalculationProcess primal_y_plus_process(
        primal_model, empty_y_plus_parameters);

    Parameters empty_nut_parameters = Parameters(R"({
        "model_part_name" : "test"
    })");
    RansNutKEpsilonHighReSensitivitiesProcess nut_sensitivities_process(
        adjoint_model, empty_nut_parameters);
    RansNutKEpsilonHighReCalculationProcess adjoint_nut_process(
        adjoint_model, empty_nut_parameters);
    RansNutKEpsilonHighReCalculationProcess primal_nut_process(primal_model, empty_nut_parameters);

    auto perturb_variable = [](NodeType& rNode) -> double& {
        return rNode.FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE_2);
    };

    auto calculate_sensitivity_matrix = [](Matrix& rOutput, Element& rElement,
                                           ProcessInfo& rProcessInfo) {
        rElement.CalculateSecondDerivativesLHS(rOutput, rProcessInfo);
        const double bossak_alpha = rProcessInfo[BOSSAK_ALPHA];
        noalias(rOutput) = rOutput * (1.0 - bossak_alpha);
    };

    RansModellingApplicationTestUtilities::RunResidualScalarSensitivityTest(
        r_primal_model_part, r_adjoint_model_part, primal_y_plus_process, primal_nut_process,
        adjoint_y_plus_process, adjoint_nut_process, y_plus_sensitivities_process,
        nut_sensitivities_process, RansEvmKEpsilonModel::UpdateVariablesInModelPart,
        calculate_sensitivity_matrix, perturb_variable, 1e-5, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmEpsilonAdjoint2D3N_CalculateSensitivityMatrix,
                          RANSModellingApplicationElementInterfaces)
{
    Model primal_model;
    ModelPart& r_primal_model_part = primal_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonElementTestModelPart(
        r_primal_model_part, "RansEvmEpsilon2D3N");

    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonElementTestModelPart(
        r_adjoint_model_part, "RansEvmEpsilonAdjoint2D3N");

    Parameters empty_y_plus_parameters = Parameters(R"({
        "model_part_name" : "test"
    })");
    RansLogarithmicYPlusVelocitySensitivitiesProcess y_plus_sensitivities_process(
        adjoint_model, empty_y_plus_parameters);
    RansLogarithmicYPlusCalculationProcess adjoint_y_plus_process(
        adjoint_model, empty_y_plus_parameters);
    RansLogarithmicYPlusCalculationProcess primal_y_plus_process(
        primal_model, empty_y_plus_parameters);

    Parameters empty_nut_parameters = Parameters(R"({
        "model_part_name" : "test"
    })");
    RansNutKEpsilonHighReSensitivitiesProcess nut_sensitivities_process(
        adjoint_model, empty_nut_parameters);
    RansNutKEpsilonHighReCalculationProcess adjoint_nut_process(
        adjoint_model, empty_nut_parameters);
    RansNutKEpsilonHighReCalculationProcess primal_nut_process(primal_model, empty_nut_parameters);

    auto perturb_variable = [](NodeType& rNode, const int Dim) -> double& {
        array_1d<double, 3>& r_coordinates = rNode.Coordinates();
        return r_coordinates[Dim];
    };

    auto calculate_sensitivity_matrix = [](Matrix& rOutput, Element& rElement,
                                           ProcessInfo& rProcessInfo) {
        rElement.CalculateSensitivityMatrix(SHAPE_SENSITIVITY, rOutput, rProcessInfo);
    };

    RansModellingApplicationTestUtilities::RunResidualVectorSensitivityTest(
        r_primal_model_part, r_adjoint_model_part, primal_y_plus_process, primal_nut_process,
        adjoint_y_plus_process, adjoint_nut_process, y_plus_sensitivities_process,
        nut_sensitivities_process, RansEvmKEpsilonModel::UpdateVariablesInModelPart,
        calculate_sensitivity_matrix, perturb_variable, 1e-7, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmEpsilonAdjoint2D3N_Calculate_RANS_VELOCITY_PRESSURE_PARTIAL_DERIVATIVE,
                          RANSModellingApplicationElementInterfaces)
{
    Model primal_model;
    ModelPart& r_primal_model_part = primal_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonElementTestModelPart(
        r_primal_model_part, "RansEvmEpsilon2D3N");

    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonElementTestModelPart(
        r_adjoint_model_part, "RansEvmEpsilonAdjoint2D3N");

    Parameters empty_y_plus_parameters = Parameters(R"({
        "model_part_name" : "test"
    })");
    RansLogarithmicYPlusVelocitySensitivitiesProcess y_plus_sensitivities_process(
        adjoint_model, empty_y_plus_parameters);
    RansLogarithmicYPlusCalculationProcess adjoint_y_plus_process(
        adjoint_model, empty_y_plus_parameters);
    RansLogarithmicYPlusCalculationProcess primal_y_plus_process(
        primal_model, empty_y_plus_parameters);

    Parameters empty_nut_parameters = Parameters(R"({
        "model_part_name" : "test"
    })");
    RansNutKEpsilonHighReSensitivitiesProcess nut_sensitivities_process(
        adjoint_model, empty_nut_parameters);
    RansNutKEpsilonHighReCalculationProcess adjoint_nut_process(
        adjoint_model, empty_nut_parameters);
    RansNutKEpsilonHighReCalculationProcess primal_nut_process(primal_model, empty_nut_parameters);

    auto perturb_variable = [](NodeType& rNode, const int Dim) -> double& {
        array_1d<double, 3>& r_velocity = rNode.FastGetSolutionStepValue(VELOCITY);
        return r_velocity[Dim];
    };

    auto calculate_sensitivity_matrix = [](Matrix& rOutput, Element& rElement,
                                           ProcessInfo& rProcessInfo) {
        rElement.Calculate(RANS_VELOCITY_PRESSURE_PARTIAL_DERIVATIVE, rOutput, rProcessInfo);
    };

    RansModellingApplicationTestUtilities::RunResidualVectorSensitivityTest(
        r_primal_model_part, r_adjoint_model_part, primal_y_plus_process, primal_nut_process,
        adjoint_y_plus_process, adjoint_nut_process, y_plus_sensitivities_process,
        nut_sensitivities_process, RansEvmKEpsilonModel::UpdateVariablesInModelPart,
        calculate_sensitivity_matrix, perturb_variable, 1e-6, 1e-5);
}

} // namespace Testing
} // namespace Kratos