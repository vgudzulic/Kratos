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

#if !defined(KRATOS_COUPLED_STRATEGY_ITEM_H)
#define KRATOS_COUPLED_STRATEGY_ITEM_H

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
class CoupledStrategyItem
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of CoupledStrategyItem
    KRATOS_CLASS_POINTER_DEFINITION(CoupledStrategyItem);

    using BaseType = SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>;

    ///@}
    ///@name Life Cycle
    ///@{

    CoupledStrategyItem(typename BaseType::Pointer pStrategy,
                        const std::string& rStrategyName,
                        const int EchoLevel)
        : mpStrategy(pStrategy),
          mStrategyName(rStrategyName),
          mEchoLevel(EchoLevel)
    {
    }

    /// Destructor.
    ~CoupledStrategyItem() = default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void AddAuxiliaryProcess(Process::Pointer pAuxiliaryProcess)
    {
        mAuxiliaryProcessList.push_back(pAuxiliaryProcess);
        KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
            << "Added " << pAuxiliaryProcess->Info() << " auxiliary process\n";
    }

    std::string GetName() const
    {
        return mStrategyName;
    }

    std::vector<Process::Pointer>& GetAuxiliaryProcessList()
    {
        return mAuxiliaryProcessList;
    }

    typename BaseType::Pointer pGetStrategy()
    {
        return mpStrategy;
    }

    BaseType& GetStrategy()
    {
        return *mpStrategy;
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
    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "CoupledStrategyItem ( " << mStrategyName << " )";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const
    {
    }

    std::string GetStrategyInfo() const
    {
        std::stringstream buffer;

        const ModelPart& r_model_part = mpStrategy->GetModelPart();

        buffer << "  Strategy name        : " << mStrategyName << "\n";
        buffer << "  Model part name      : " << r_model_part.Name() << "\n";
        buffer << "  Element type         : " << r_model_part.ElementsBegin()->Info() << "\n";
        buffer << "  Condition type       : " << r_model_part.ConditionsBegin()->Info()
               << "\n";

        buffer << "  Update process list:\n";
        if (mAuxiliaryProcessList.size() == 0)
        {
            buffer << "      No update processes.\n";
        }
        else
        {
            for (const Process::Pointer& r_process_pointer : mAuxiliaryProcessList)
            {
                buffer << "      Process name   : " << r_process_pointer->Info() << "\n";
            }
        }

        return buffer.str();
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

    typename BaseType::Pointer mpStrategy;
    std::string mStrategyName;
    int mEchoLevel;

    std::vector<Process::Pointer> mAuxiliaryProcessList;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

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
    CoupledStrategyItem& operator=(CoupledStrategyItem const& rOther)
    {
    }

    /// Copy constructor.
    CoupledStrategyItem(CoupledStrategyItem const& rOther)
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

#endif // KRATOS_COUPLED_STRATEGY_ITEM_H
