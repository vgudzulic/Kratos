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

#if !defined(KRATOS_DENSE_LU_DIRECT_SOLVER_H_INCLUDED)
#define KRATOS_DENSE_LU_DIRECT_SOLVER_H_INCLUDED

// External includes

// Project includes
#include "includes/define.h"
#include "linear_solvers/direct_solver.h"
#include "utilities/math_utils.h"
#include "factories/dense_linear_solver_factory.h"

namespace Kratos {

template <
    class TSparseSpaceType,
    class TDenseSpaceType,
    class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType>>
class DenseLUDirectSolver
    : public DirectSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>
{
private:
    DenseLUDirectSolver &operator=(const DenseLUDirectSolver &Other);

    DenseLUDirectSolver(const DenseLUDirectSolver &Other);

public:
    KRATOS_CLASS_POINTER_DEFINITION(DenseLUDirectSolver);

    typedef DirectSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType MatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TSparseSpaceType::DataType DataType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    DenseLUDirectSolver() {}

    DenseLUDirectSolver(Parameters settings) : BaseType(settings)
    {
    }

    ~DenseLUDirectSolver() override {}

    /**
     * @brief This function is designed to be called every time the coefficients change in the system that is, normally at the beginning of each solve.
     * @details For example if we are implementing a direct solver, this is the place to do the factorization
     * so that then the backward substitution can be performed effectively more than once
     * @param rA System matrix
     * @param rX Solution vector
     * @param rB Right hand side vector
     */
    void InitializeSolutionStep(MatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        // NOTE: DOES NOTHING
    }

    /**
     * @brief This function actually performs the solution work, eventually taking advantage of what was done before in the
     * @details Initialize and InitializeSolutionStep functions.
     * @param rA System matrix
     * @param rX Solution vector
     * @param rB Right hand side vector
     */
    void PerformSolutionStep(MatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        double det;
        MatrixType invA;
        MathUtils<double>::InvertMatrix(rA, invA, det);
        if (rX.size() != rB.size()) {
            rX.resize(rB.size());
        }
        noalias(rX) = prod(invA, rB);
    }

    /**
     * @brief Solves the linear system Ax=b
     * @param rA System matrix
     * @param rX Solution vector
     * @param rB Right hand side vector
     * @return true if solution found, otherwise false
     */
    bool Solve(MatrixType &rA, VectorType &rX, VectorType &rB) override
    {
        InitializeSolutionStep(rA, rX, rB);
        PerformSolutionStep(rA, rX, rB);

        return true;
    }

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "DenseLUDirectSolver";
        return  buffer.str();
    }

    /**
     * Print information about this object.
     */
    void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << Info();
    }

    /**
     * Print object's data.
     */
    void PrintData(std::ostream &rOStream) const override
    {
    }

    static DenseLinearSolverFactory<TSparseSpaceType, TDenseSpaceType, DenseLUDirectSolver> Factory()
    {
        return DenseLinearSolverFactory<TSparseSpaceType, TDenseSpaceType, DenseLUDirectSolver>();
    }
}; // class DenseLUDirectSolver

/**
 * input stream function
 */
template<
    class TSparseSpaceType,
    class TDenseSpaceType,
    class TReordererType>
inline std::istream &operator>>(
    std::istream &rIStream,
    DenseLUDirectSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> &rThis)
{
    return rIStream;
}

/**
 * output stream function
 */
template<
    class TSparseSpaceType,
    class TDenseSpaceType,
    class TReordererType
    >
inline std::ostream &operator<<(
    std::ostream &rOStream,
    const DenseLUDirectSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> &rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos

#endif // defined(KRATOS_DENSE_LU_DIRECT_SOLVER_H_INCLUDED)
