//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes
#include <cmath>
#include <limits>

// External includes

// Project includes
#include "includes/define.h"
#include "utilities/variable_utils.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

// Application includes
#include "rans_application_variables.h"

// Include base h
#include "rans_variable_utilities.h"

namespace Kratos
{
namespace RansVariableUtilities
{
std::tuple<unsigned int, unsigned int> ClipScalarVariable(
    const double MinimumValue,
    const double MaximumValue,
    const Variable<double>& rVariable,
    ModelPart& rModelPart)
{
    KRATOS_TRY

    auto& r_communicator = rModelPart.GetCommunicator();
    auto& r_nodes = r_communicator.LocalMesh().Nodes();

    unsigned int number_of_nodes_below_minimum, number_of_nodes_above_maximum;
    std::tie(number_of_nodes_below_minimum, number_of_nodes_above_maximum) =
        BlockPartition<ModelPart::NodesContainerType>(r_nodes)
            .for_each<CombinedReduction<SumReduction<unsigned int>, SumReduction<unsigned int>>>(
                [&](ModelPart::NodeType& rNode) -> std::tuple<unsigned int, unsigned int> {
                    double& r_value = rNode.FastGetSolutionStepValue(rVariable);

                    if (r_value < MinimumValue) {
                        r_value = MinimumValue;
                        return std::tuple<unsigned int, unsigned int>(1, 0);
                    } else if (r_value > MaximumValue) {
                        r_value = MaximumValue;
                        return std::tuple<unsigned int, unsigned int>(0, 1);
                    }

                    return std::tuple<unsigned int, unsigned int>(0, 0);
                });

    r_communicator.SynchronizeVariable(rVariable);

    // Stores followings
    // index - 0 : number_of_nodes_below_minimum
    // index - 1 : number_of_nodes_above_maximum
    std::vector<unsigned int> nodes_count = {number_of_nodes_below_minimum,
                                             number_of_nodes_above_maximum};
    const std::vector<unsigned int>& total_nodes_count =
        r_communicator.GetDataCommunicator().SumAll(nodes_count);

    return std::tuple<unsigned int, unsigned int>(total_nodes_count[0], total_nodes_count[1]);

    KRATOS_CATCH("")
}

double GetMinimumScalarValue(
    const ModelPart& rModelPart,
    const Variable<double>& rVariable)
{
    KRATOS_TRY

    const auto& r_communicator = rModelPart.GetCommunicator();
    const auto& r_nodes = r_communicator.LocalMesh().Nodes();

    const int number_of_nodes = r_nodes.size();
    const double min_value =
        IndexPartition<int>(number_of_nodes).for_each<MinReduction<double>>([&](const int i) -> double {
            return (r_nodes.begin() + i)->FastGetSolutionStepValue(rVariable);
        });

    return r_communicator.GetDataCommunicator().MinAll(min_value);

    KRATOS_CATCH("");
}

double GetMaximumScalarValue(
    const ModelPart& rModelPart,
    const Variable<double>& rVariable)
{
    KRATOS_TRY

    const auto& r_communicator = rModelPart.GetCommunicator();
    const auto& r_nodes = r_communicator.LocalMesh().Nodes();

    const int number_of_nodes = r_nodes.size();
    const double max_value =
        IndexPartition<int>(number_of_nodes).for_each<MaxReduction<double>>([&](const int i) -> double {
            return (r_nodes.begin() + i)->FastGetSolutionStepValue(rVariable);
        });

    return r_communicator.GetDataCommunicator().MaxAll(max_value);

    KRATOS_CATCH("");
}

void GetNodalVariablesVector(
    Vector& rValues,
    const ModelPart::NodesContainerType& rNodes,
    const Variable<double>& rVariable)
{
    const int number_of_nodes = rNodes.size();

    if (static_cast<int>(rValues.size()) != number_of_nodes) {
        rValues.resize(number_of_nodes);
    }

    IndexPartition<int>(number_of_nodes).for_each([&](const int i_node) {
        rValues[i_node] = (rNodes.begin() + i_node)->FastGetSolutionStepValue(rVariable);
    });
}

void GetNodalArray(
    Vector& rNodalValues,
    const Element& rElement,
    const Variable<double>& rVariable)
{
    const auto& r_geometry = rElement.GetGeometry();
    std::size_t number_of_nodes = r_geometry.PointsNumber();

    if (rNodalValues.size() != number_of_nodes) {
        rNodalValues.resize(number_of_nodes);
    }

    for (std::size_t i_node = 0; i_node < number_of_nodes; ++i_node) {
        rNodalValues[i_node] = r_geometry[i_node].FastGetSolutionStepValue(rVariable);
    }
}

void SetNodalVariables(
    ModelPart::NodesContainerType& rNodes,
    const Vector& rValues,
    const Variable<double>& rVariable)
{
    const int number_of_nodes = rNodes.size();

    KRATOS_ERROR_IF(static_cast<int>(rValues.size()) != number_of_nodes)
        << "rValues vector size mismatch with rNodes size in "
           "SetNodalVariables. [ rValues.size = "
        << rValues.size() << ", rNodes.size = " << rNodes.size() << " ]\n";

    IndexPartition<int>(number_of_nodes).for_each([&](const int i_node) {
        (rNodes.begin() + i_node)->FastGetSolutionStepValue(rVariable) = rValues[i_node];
    });
}

void CopyNodalSolutionStepVariablesList(
    ModelPart& rOriginModelPart,
    ModelPart& rDestinationModelPart)
{
    KRATOS_TRY

    rDestinationModelPart.GetNodalSolutionStepVariablesList() =
        rOriginModelPart.GetNodalSolutionStepVariablesList();

    KRATOS_CATCH("");
}

template <typename TDataType>
void AssignConditionVariableValuesToNodes(ModelPart& rModelPart,
                                          const Variable<TDataType>& rVariable,
                                          const Flags& rFlag,
                                          const bool FlagValue)
{
    auto& r_nodes = rModelPart.Nodes();
    VariableUtils().SetHistoricalVariableToZero(rVariable, r_nodes);

    const int number_of_conditions = rModelPart.NumberOfConditions();
#pragma omp parallel for
    for (int i_cond = 0; i_cond < number_of_conditions; ++i_cond)
    {
        auto& r_cond = *(rModelPart.ConditionsBegin() + i_cond);
        if (r_cond.Is(rFlag) == FlagValue)
        {
            const int number_of_nodes = r_cond.GetGeometry().PointsNumber();
            const auto& r_normal = r_cond.GetValue(rVariable);
            for (int i_node = 0; i_node < number_of_nodes; ++i_node)
            {
                auto& r_node = r_cond.GetGeometry()[i_node];
                r_node.SetLock();
                r_node.FastGetSolutionStepValue(rVariable) +=
                    r_normal * (1.0 / static_cast<double>(number_of_nodes));
                r_node.UnSetLock();
            }
        }
    }

    rModelPart.GetCommunicator().AssembleCurrentData(rVariable);
}

void AddAnalysisStep(ModelPart& rModelPart, const std::string& rStepName)
{
    auto& r_process_info = rModelPart.GetProcessInfo();
    if (!r_process_info.Has(ANALYSIS_STEPS))
    {
        r_process_info.SetValue(ANALYSIS_STEPS, std::vector<std::string>());
    }
    r_process_info[ANALYSIS_STEPS].push_back(rStepName);
}

bool IsAnalysisStepCompleted(const ModelPart& rModelPart, const std::string& rStepName)
{
    const auto& r_process_info = rModelPart.GetProcessInfo();
    if (r_process_info.Has(ANALYSIS_STEPS))
    {
        const std::vector<std::string>& r_steps = r_process_info[ANALYSIS_STEPS];
        return (std::find(r_steps.begin(), r_steps.end(), rStepName) != r_steps.end());
    }
    else
    {
        return false;
    }
}

void AssignBoundaryFlagsToGeometries(ModelPart& rModelPart)
{
    const int number_of_conditions = rModelPart.NumberOfConditions();
#pragma omp parallel for
    for (int i_cond = 0; i_cond < number_of_conditions; ++i_cond)
    {
        Condition& r_condition = *(rModelPart.ConditionsBegin() + i_cond);
        r_condition.SetValue(RANS_IS_INLET, r_condition.Is(INLET));
        r_condition.SetValue(RANS_IS_OUTLET, r_condition.Is(OUTLET));
        r_condition.SetValue(RANS_IS_STRUCTURE, r_condition.Is(STRUCTURE));
    }
}

void FixFlaggedDofs(ModelPart& rModelPart,
                    const Variable<double>& rFixingVariable,
                    const Flags& rFlag,
                    const bool CheckValue)
{
    KRATOS_TRY

    const int number_of_nodes = rModelPart.NumberOfNodes();
#pragma omp parallel for
    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        auto& r_node = *(rModelPart.NodesBegin() + i_node);
        if (r_node.Is(rFlag) == CheckValue)
        {
            r_node.Fix(rFixingVariable);
        }
    }

    KRATOS_CATCH("");
}

void CalculateMagnitudeSquareForNodal3DVariable(ModelPart& rModelPart,
                                                const Variable<array_1d<double, 3>>& r3DVariable,
                                                const Variable<double>& rOutputVariable)
{
    const int number_of_nodes = rModelPart.NumberOfNodes();
#pragma omp parallel for
    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        auto& r_node = *(rModelPart.NodesBegin() + i_node);
        const double magnitude = norm_2(r_node.FastGetSolutionStepValue(r3DVariable));
        r_node.FastGetSolutionStepValue(rOutputVariable) = std::pow(magnitude, 2);
    }
}

template <>
KRATOS_API(RANS_APPLICATION)
double GetVariableValueNorm(const double& rValue)
{
    return rValue;
}

template <>
KRATOS_API(RANS_APPLICATION)
double GetVariableValueNorm(const array_1d<double, 3>& rValue)
{
    return norm_2(rValue);
}

template <typename TDataType>
std::tuple<double, double> CalculateTransientVariableConvergence(const ModelPart& rModelPart,
                                                                 const Variable<TDataType>& rVariable)
{
    KRATOS_TRY

    const auto& r_communicator = rModelPart.GetCommunicator();
    const auto& r_nodes = r_communicator.LocalMesh().Nodes();
    const int number_of_nodes = r_nodes.size();

    KRATOS_ERROR_IF(rModelPart.GetBufferSize() < 2)
        << rModelPart.Name() << " buffer size is "
        << rModelPart.GetBufferSize() << ". Buffer size of 2 or greater is required to calculate transient variable convergence for "
        << rVariable.Name() << ".\n";

    double dx = 0.0;
    double number_of_dofs = 0.0;
    double solution = 0.0;
#pragma omp parallel for reduction(+ : dx, number_of_dofs, solution)
    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        const auto& r_node = *(r_nodes.begin() + i_node);
        const auto& r_old_value = r_node.FastGetSolutionStepValue(rVariable, 1);
        const auto& r_new_value = r_node.FastGetSolutionStepValue(rVariable);
        dx += std::pow(GetVariableValueNorm<TDataType>(r_new_value - r_old_value), 2);
        solution += std::pow(GetVariableValueNorm<TDataType>(r_new_value), 2);
        number_of_dofs += ((r_node.HasDofFor(rVariable) && !r_node.IsFixed(rVariable)) ||
                           !r_node.HasDofFor(rVariable));
    }

    // to improve mpi communication performance
    const std::vector<double> process_values = {dx, solution, number_of_dofs};
    const std::vector<double>& global_values =
        r_communicator.GetDataCommunicator().SumAll(process_values);

    dx = std::sqrt(global_values[0]);
    solution = std::sqrt(global_values[1]);
    number_of_dofs = std::max(1.0, global_values[2]);
    solution = (solution > 0.0 ? solution : 1.0);

    return std::make_tuple<double, double>(std::forward<double>(dx / solution),
                                           std::forward<double>(dx / number_of_dofs));

    KRATOS_CATCH("");
}

// template instantiations
template KRATOS_API(RANS_APPLICATION) void AssignConditionVariableValuesToNodes<double>(
    ModelPart&, const Variable<double>&, const Flags&, const bool);

template KRATOS_API(RANS_APPLICATION) void AssignConditionVariableValuesToNodes<array_1d<double, 3>>(
    ModelPart&, const Variable<array_1d<double, 3>>&, const Flags&, const bool);

template KRATOS_API(RANS_APPLICATION) std::tuple<double, double> CalculateTransientVariableConvergence<double>(
    const ModelPart&, const Variable<double>&);

template KRATOS_API(RANS_APPLICATION)
    std::tuple<double, double> CalculateTransientVariableConvergence<array_1d<double, 3>>(
        const ModelPart&, const Variable<array_1d<double, 3>>&);

} // namespace RansVariableUtilities

} /* namespace Kratos.*/
