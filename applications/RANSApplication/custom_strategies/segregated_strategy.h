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
#include "rans_application_variables.h"

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

    void AddStrategy(typename BaseType::Pointer pStrategy, const std::string& rStrategyName)
    {
        mSolvingStrategiesList.push_back(pStrategy);
        mSolvingStrategyNamesList.push_back(rStrategyName);
        mIsStrategyConverged.push_back(false);
    }

    void AddAuxiliaryProcess(Process::Pointer pAuxiliaryProcess, const std::string& rVariableName)
    {
        KRATOS_TRY

        const int number_of_processes = mSolvingStrategiesList.size();
        int current_strategy_index = -1;
        for (int i = 0; i < number_of_processes; ++i)
        {
            if (mSolvingStrategyNamesList[i] == rVariableName)
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
        int iteration_format_length =
            static_cast<int>(std::log10(this->mMaxIterations)) + 1;

        const ProcessInfo& r_current_process_info = this->GetModelPart().GetProcessInfo();
        const int iteration = r_current_process_info[COUPLING_ITERATION];

        bool is_converged = true;
        const int number_of_solving_strategies = this->mSolvingStrategiesList.size();

        std::stringstream conv_msg;
        conv_msg << "[Itr.#" << std::setw(iteration_format_length) << iteration
                 << "/" << this->mMaxIterations << "] CONVERGENCE CHECK: \n";

        for (int i = 0; i < number_of_solving_strategies; ++i)
        {
            auto p_solving_strategy = this->mSolvingStrategiesList[i];
            const auto& strategy_name = this->mSolvingStrategyNamesList[i];

            if (!p_solving_strategy->IsConverged())
            {
                is_converged = false;
                mIsStrategyConverged[i] = false;
                conv_msg << "      " << strategy_name << " : Not Converged\n";
            }
            else
            {
                mIsStrategyConverged[i] = true;
                conv_msg << "      " << strategy_name << " : Converged\n";
            }
        }

        KRATOS_INFO_IF(this->Info(), this->GetEchoLevel() > 1) << conv_msg.str();

        return is_converged;
    }

    bool SolveSolutionStep() override
    {
        KRATOS_TRY

        ModelPart& r_model_part = this->GetModelPart();
        ProcessInfo& r_current_process_info = r_model_part.GetProcessInfo();

        const int number_of_solving_strategies = this->mSolvingStrategiesList.size();

        bool is_converged = false;
        int& iteration = r_current_process_info[COUPLING_ITERATION];
        iteration = 1;

        while (!is_converged && iteration <= this->mMaxIterations)
        {
            for (int i = 0; i < number_of_solving_strategies; ++i)
            {
                auto p_solving_strategy = this->mSolvingStrategiesList[i];
                auto variable_name = this->mSolvingStrategyNamesList[i];

                if (!mIsStrategyConverged[i])
                {
                    p_solving_strategy->SolveSolutionStep();
                    ExecuteStrategyUpdateProcesses(i);

                    const unsigned int iterations =
                        r_current_process_info[NL_ITERATION_NUMBER];
                    KRATOS_INFO_IF(this->Info(), this->GetEchoLevel() > 1)
                        << "Solving " << variable_name << " used " << iterations
                        << " non-linear iterations.\n";
                }
            }

            is_converged = this->IsConverged();
            ++iteration;
        }

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

        for (auto p_solving_strategy : this->mSolvingStrategiesList)
            p_solving_strategy->InitializeSolutionStep();

        KRATOS_CATCH("");
    }

    void Predict() override
    {
        for (auto p_solving_strategy : this->mSolvingStrategiesList)
            p_solving_strategy->Predict();
    }

    void FinalizeSolutionStep() override
    {
        KRATOS_TRY

        for (auto p_solving_strategy : this->mSolvingStrategiesList)
            p_solving_strategy->FinalizeSolutionStep();

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

    int mMaxIterations;
    int iteration;

    std::vector<typename BaseType::Pointer> mSolvingStrategiesList;
    std::vector<std::string> mSolvingStrategyNamesList;
    std::vector<Process::Pointer> mAuxiliaryProcessList;
    std::vector<int> mAuxiliaryProcessStrategyId;
    std::vector<bool> mIsStrategyConverged;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

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
            << mSolvingStrategyNamesList[StrategyId] << ".\n";

        KRATOS_CATCH("");
    }

    std::string PrintStrategyInfo()
    {
        std::stringstream buffer;
        buffer << "SegregatedStrategy with followings:\n";
        buffer << "    Strategies list:\n";

        const int number_of_strategies = mSolvingStrategyNamesList.size();
        for (int i = 0; i < number_of_strategies; ++i)
        {
            buffer << "        Strategy id    : " << i << "\n";
            buffer << "            Strategy name     : "
                   << mSolvingStrategyNamesList[i] << "\n";
            buffer << "            Model part name   : "
                   << mSolvingStrategiesList[i]->GetModelPart().Name() << "\n";
            buffer << "            Element type      : "
                   << mSolvingStrategiesList[i]->GetModelPart().ElementsBegin()->Info()
                   << "\n";
            buffer
                << "            Condition type    : "
                << mSolvingStrategiesList[i]->GetModelPart().ConditionsBegin()->Info()
                << "\n";
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
