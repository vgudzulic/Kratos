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

#if !defined(KRATOS_RANS_SEGREGATED_STRATEGY_H)
#define KRATOS_RANS_SEGREGATED_STRATEGY_H

// System includes
#include <string>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "solving_strategies/strategies/solving_strategy.h"

// Application includes
#include "custom_utilities/rans_variable_utilities.h"

namespace Kratos
{
///@addtogroup RANSApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}

///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
class SegregatedStrategy : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of SegregatedStrategy
    KRATOS_CLASS_POINTER_DEFINITION(SegregatedStrategy);

    using BaseType = SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>;

    ///@}
    ///@name Life Cycle
    ///@{

    SegregatedStrategy(ModelPart& rModelPart, const int MaxIterations = 10)
        : BaseType(rModelPart, false), mMaxIterations(MaxIterations)
    {
        KRATOS_TRY;

        KRATOS_CATCH("");
    }

    /// Destructor.
    ~SegregatedStrategy() override = default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void AddStrategy(typename BaseType::Pointer pStrategy,
                     const std::string& rVariableName,
                     const double RelativeTolerance = 1e-3,
                     const double AbsoluteTolerance = 1e-5)
    {
        mSolvingStrategiesList.push_back(pStrategy);
        mSolvingVariableNamesList.push_back(rVariableName);
        mSolvingVariableRelativeToleranceList.push_back(RelativeTolerance);
        mSolvingVariableAbsoluteToleranceList.push_back(AbsoluteTolerance);
    }

    void AddAuxiliaryProcess(Process::Pointer pAuxiliaryProcess, const std::string& rVariableName)
    {
        KRATOS_TRY

        const int number_of_processes = mSolvingStrategiesList.size();
        int current_strategy_index = -1;
        for (int i = 0; i < number_of_processes; ++i)
        {
            if (mSolvingVariableNamesList[i] == rVariableName)
            {
                current_strategy_index = i;
                break;
            }
        }

        KRATOS_ERROR_IF(current_strategy_index == -1)
            << "A strategy for " << rVariableName << " not found.\n";

        mAuxiliaryProcessList.push_back(pAuxiliaryProcess);
        mAuxiliaryProcessStrategyId.push_back(current_strategy_index);

        KRATOS_CATCH("");
    }

    void Initialize() override
    {
        KRATOS_TRY

        for (auto process : mAuxiliaryProcessList)
            process->ExecuteInitialize();

        for (auto strategy : mSolvingStrategiesList)
            strategy->Initialize();

        KRATOS_INFO_IF(this->Info(), this->GetEchoLevel() > 0)
            << this->PrintStrategyInfo();

        KRATOS_CATCH("");
    }

    int Check() override
    {
        KRATOS_TRY;

        for (auto variable_name : mSolvingVariableNamesList)
        {
            KRATOS_ERROR_IF((!KratosComponents<Variable<double>>::Has(variable_name) &&
                             !KratosComponents<Variable<array_1d<double, 3>>>::Has(variable_name)))
                << variable_name << " not found in either double or 3d variables list.\n";
        }

        for (auto process : mAuxiliaryProcessList)
            process->Check();

        for (auto strategy : mSolvingStrategiesList)
            strategy->Check();

        KRATOS_ERROR_IF(mSolvingStrategiesList.size() == 0)
            << "No strategies are found for SegregatedStrategy.";

        return 0;

        KRATOS_CATCH("");
    }

    bool IsConverged() override
    {
        return true;
    }

    bool SolveSolutionStep() override
    {
        KRATOS_TRY

        ModelPart& r_model_part = this->GetModelPart();
        Communicator& r_communicator = r_model_part.GetCommunicator();
        ProcessInfo& r_current_process_info = r_model_part.GetProcessInfo();
        ModelPart::NodesContainerType& r_nodes = r_communicator.LocalMesh().Nodes();

        const int number_of_solving_strategies = this->mSolvingStrategiesList.size();

        for (auto p_solving_strategy : this->mSolvingStrategiesList)
        {
            p_solving_strategy->InitializeSolutionStep();
            p_solving_strategy->Predict();
        }

        int iteration_format_length =
            static_cast<int>(std::log10(this->mMaxIterations)) + 1;

        Vector old_values;
        Vector new_values;
        Vector delta_values;

        bool is_converged = false;
        int iteration = 1;

        std::vector<int> is_strategy_converged;
        is_strategy_converged.resize(number_of_solving_strategies);

        while (!is_converged && iteration <= this->mMaxIterations)
        {
            for (int i = 0; i < number_of_solving_strategies; ++i)
            {
                auto p_solving_strategy = this->mSolvingStrategiesList[i];
                auto variable_name = this->mSolvingVariableNamesList[i];

                GetNodalVariablesVector(old_values, r_nodes, variable_name);

                p_solving_strategy->SolveSolutionStep();
                ExecuteStrategyUpdateProcesses(i);

                GetNodalVariablesVector(new_values, r_nodes, variable_name);

                if (delta_values.size() < new_values.size())
                    delta_values.resize(new_values.size());
                noalias(delta_values) = new_values - old_values;

                // This vector stores norms of the residual
                // index - 0 : increase_norm
                // index - 1 : solution_norm
                // index - 3 : number of nodes
                std::vector<double> residual_norms(3);
                residual_norms[0] = std::pow(norm_2(delta_values), 2);
                residual_norms[1] = std::pow(norm_2(new_values), 2);
                residual_norms[2] = static_cast<double>(r_nodes.size());
                const std::vector<double>& total_residual_norms =
                    r_communicator.GetDataCommunicator().SumAll(residual_norms);

                double convergence_relative =
                    total_residual_norms[0] /
                    (total_residual_norms[1] <= std::numeric_limits<double>::epsilon()
                         ? 1.0
                         : total_residual_norms[1]);
                double convergence_absolute =
                    std::sqrt(total_residual_norms[0]) / total_residual_norms[2];

                is_strategy_converged[i] =
                    (convergence_relative < this->mSolvingVariableRelativeToleranceList[i] ||
                     convergence_absolute < this->mSolvingVariableAbsoluteToleranceList[i]);

                if (this->GetEchoLevel() > 1)
                {
                    const unsigned int iterations =
                        r_current_process_info[NL_ITERATION_NUMBER];
                    KRATOS_INFO(this->Info())
                        << "Solving " << variable_name << " used " << iterations
                        << " non-linear iterations.\n";

                    std::stringstream conv_check_msg;
                    conv_check_msg
                        << "[Itr.#" << std::setw(iteration_format_length)
                        << iteration << "/" << this->mMaxIterations
                        << "] CONVERGENCE CHECK: " << variable_name
                        << " ratio = " << std::setprecision(3) << std::scientific
                        << convergence_relative << "; exp. ratio = "
                        << mSolvingVariableRelativeToleranceList[i]
                        << "; abs = " << convergence_absolute << "; exp.abs = "
                        << mSolvingVariableAbsoluteToleranceList[i] << "\n";
                    KRATOS_INFO_IF(this->Info(), this->GetEchoLevel() > 2)
                        << conv_check_msg.str();

                    if (is_strategy_converged[i])
                    {
                        std::stringstream conv_msg;
                        conv_msg << "[Itr.#" << std::setw(iteration_format_length)
                                 << iteration << "/" << this->mMaxIterations
                                 << "] CONVERGENCE CHECK: " << variable_name
                                 << " *** CONVERGENCE IS ACHIEVED ***\n";
                        KRATOS_INFO(this->Info()) << conv_msg.str();
                    }
                }
            }

            is_converged = true;
            for (int i = 0; i < number_of_solving_strategies; ++i)
            {
                is_converged = std::min((int)(is_converged), is_strategy_converged[i]);
            }
        }

        for (auto p_solving_strategy : this->mSolvingStrategiesList)
            p_solving_strategy->FinalizeSolutionStep();

        KRATOS_INFO_IF(this->Info(), !is_converged && this->GetEchoLevel() > 0)
            << "\n-------------------------------------------------------"
            << "\n    INFO: Max coupling iterations reached.             "
            << "\n          Please increase coupling max_iterations      "
            << "\n          or decrease coupling                         "
            << "\n          relative_tolerance/absolute tolerance        "
            << "\n-------------------------------------------------------"
            << "\n";

        return is_converged;

        KRATOS_CATCH("");
    }

    void InitializeSolutionStep() override
    {
        KRATOS_TRY

        for (auto process : mAuxiliaryProcessList)
            process->ExecuteInitializeSolutionStep();

        KRATOS_CATCH("");
    }

    void FinalizeSolutionStep() override
    {
        KRATOS_TRY

        for (auto process : mAuxiliaryProcessList)
            process->ExecuteFinalizeSolutionStep();

        KRATOS_CATCH("");
    }

    void Clear() override
    {
        for (auto strategy : mSolvingStrategiesList)
            strategy->Clear();
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "SegregatedStrategy";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    ///@name Protected Life Cycle
    ///@{

    ///@}
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    double mMaxIterations;

    std::vector<typename BaseType::Pointer> mSolvingStrategiesList;
    std::vector<std::string> mSolvingVariableNamesList;
    std::vector<double> mSolvingVariableRelativeToleranceList;
    std::vector<double> mSolvingVariableAbsoluteToleranceList;
    std::vector<Process::Pointer> mAuxiliaryProcessList;
    std::vector<int> mAuxiliaryProcessStrategyId;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void GetNodalVariablesVector(Vector& rOutput,
                                 const ModelPart::NodesContainerType& rNodes,
                                 const std::string& rVariableName)
    {
        KRATOS_TRY

        if (KratosComponents<Variable<double>>::Has(rVariableName))
        {
            const Variable<double>& r_variable =
                KratosComponents<Variable<double>>::Get(rVariableName);
            RansVariableUtilities::GetNodalVariablesVector(rOutput, rNodes, r_variable);
        }
        else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(rVariableName))
        {
            const Variable<array_1d<double, 3>>& r_variable =
                KratosComponents<Variable<array_1d<double, 3>>>::Get(rVariableName);
            RansVariableUtilities::GetNodalVariablesVector(rOutput, rNodes, r_variable);
        }
        else
        {
            KRATOS_ERROR << "Unsupported " << rVariableName << ".";
        }

        KRATOS_CATCH("");
    }

    void ExecuteStrategyUpdateProcesses(const int StrategyId)
    {
        KRATOS_TRY

        const int number_of_processes = mAuxiliaryProcessList.size();
        for (int i = 0; i < number_of_processes; ++i)
        {
            if (mAuxiliaryProcessStrategyId[i] == StrategyId)
            {
                mAuxiliaryProcessList[i]->Execute();
            }
        }

        KRATOS_INFO_IF(this->Info(), this->GetEchoLevel() > 1)
            << "Executed update processes for "
            << mSolvingVariableNamesList[StrategyId] << ".\n";

        KRATOS_CATCH("");
    }

    std::string PrintStrategyInfo()
    {
        std::stringstream buffer;
        buffer << "SegregatedStrategy with followings:\n";
        buffer << "    Strategies list:\n";

        const int number_of_strategies = mSolvingVariableNamesList.size();
        for (int i = 0; i < number_of_strategies; ++i)
        {
            buffer << "        Strategy id    : " << i << "\n";
            buffer << "            Variable name     : "
                   << mSolvingVariableNamesList[i] << "\n";
            buffer << "            Model part name   : "
                   << mSolvingStrategiesList[i]->GetModelPart().Name() << "\n";
            buffer << "            Element type      : "
                   << mSolvingStrategiesList[i]->GetModelPart().ElementsBegin()->Info()
                   << "\n";
            buffer
                << "            Condition type    : "
                << mSolvingStrategiesList[i]->GetModelPart().ConditionsBegin()->Info()
                << "\n";
            buffer << "            Relative tolerance: "
                   << mSolvingVariableRelativeToleranceList[i] << "\n";
            buffer << "            Absolute tolerance: "
                   << mSolvingVariableAbsoluteToleranceList[i] << "\n";
        }

        buffer << "    Update process list:\n";
        const int number_of_processes = mAuxiliaryProcessList.size();
        for (int i = 0; i < number_of_processes; ++i)
        {
            buffer << "      Process name   : " << mAuxiliaryProcessList[i]->Info() << "\n";
            buffer << "         Strategy id : " << mAuxiliaryProcessStrategyId[i] << "\n";
        }

        return buffer.str();
    }

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    SegregatedStrategy& operator=(SegregatedStrategy const& rOther)
    {
    }

    /// Copy constructor.
    SegregatedStrategy(SegregatedStrategy const& rOther)
    {
    }

    ///@}

}; /// Class FStepStrategy

///@}
///@name Type Definitions
///@{

///@}

///@} // addtogroup

} // namespace Kratos.

#endif // KRATOS_RANS_SEGREGATED_STRATEGY_H
