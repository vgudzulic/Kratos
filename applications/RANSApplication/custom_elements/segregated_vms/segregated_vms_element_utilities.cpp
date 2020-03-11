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
#include <cmath>

// Application includes

// Include base h
#include "segregated_vms_element_utilities.h"

namespace Kratos
{
namespace SegregatedVMSElementUtilities
{
template <IndexType TDim, IndexType TNumNodes>
BoundedVector<double, TNumNodes> GetConvectionOperator(const array_1d<double, 3>& rVector,
                                                       const Matrix& rShapeDerivatives)
{
    KRATOS_TRY

    KRATOS_DEBUG_ERROR_IF(TNumNodes != rShapeDerivatives.size1())
        << "Shape derivatives size1 and number of nodes mismatch [ "
        << rShapeDerivatives.size1() << " != " << TNumNodes << " ].\n";

    KRATOS_DEBUG_ERROR_IF(TDim != rShapeDerivatives.size2())
        << "Shape derivatives size2 and dimension mismatch [ "
        << rShapeDerivatives.size2() << " != " << TDim << " ].\n";

    BoundedVector<double, TNumNodes> result;
    for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) // Loop over nodes
    {
        result[i_node] = rVector[0] * rShapeDerivatives(i_node, 0);
        for (IndexType d = 1; d < TDim; ++d) // loop over components
            result[i_node] += rVector[d] * rShapeDerivatives(i_node, d);
    }

    return result;

    KRATOS_CATCH("");
}

double CalculateTauOne(const double VelocityMagnitude,
                       const double ElemSize,
                       const double Density,
                       const double Viscosity,
                       const double DynamicTau,
                       const double DeltaTime)
{
    double InvTau = Density * (DynamicTau / DeltaTime + 2.0 * VelocityMagnitude / ElemSize) +
                    4.0 * Viscosity / (ElemSize * ElemSize);
    return 1.0 / InvTau;
}

double CalculateTauTwo(const double VelocityMagnitude,
                       const double ElemSize,
                       const double Density,
                       const double Viscosity)
{
    return (Viscosity + 0.5 * Density * ElemSize * VelocityMagnitude);
}

template <>
void AddViscousTerm<2, 3>(BoundedMatrix<double, 6, 6>& rOutput,
                          const Matrix& rShapeDerivatives,
                          const double DynamicViscosity)
{
    KRATOS_TRY

    const IndexType number_of_nodes = 3;

    KRATOS_DEBUG_ERROR_IF(number_of_nodes != rShapeDerivatives.size1())
        << "Shape derivatives size1 and number of nodes mismatch [ "
        << rShapeDerivatives.size1() << " != " << number_of_nodes << " ].\n";

    KRATOS_DEBUG_ERROR_IF(rShapeDerivatives.size2() != 2)
        << "Shape derivatives size2 and dimension mismatch [ "
        << rShapeDerivatives.size2() << " != 2 ].\n";

    const double FourThirds = 4.0 / 3.0;
    const double nTwoThirds = -2.0 / 3.0;

    IndexType FirstRow(0), FirstCol(0);

    for (IndexType j = 0; j < number_of_nodes; ++j)
    {
        for (IndexType i = 0; i < number_of_nodes; ++i)
        {
            // First Row
            rOutput(FirstRow, FirstCol) +=
                DynamicViscosity *
                (FourThirds * rShapeDerivatives(i, 0) * rShapeDerivatives(j, 0) +
                 rShapeDerivatives(i, 1) * rShapeDerivatives(j, 1));
            rOutput(FirstRow, FirstCol + 1) +=
                DynamicViscosity *
                (nTwoThirds * rShapeDerivatives(i, 0) * rShapeDerivatives(j, 1) +
                 rShapeDerivatives(i, 1) * rShapeDerivatives(j, 0));

            // Second Row
            rOutput(FirstRow + 1, FirstCol) +=
                DynamicViscosity *
                (nTwoThirds * rShapeDerivatives(i, 1) * rShapeDerivatives(j, 0) +
                 rShapeDerivatives(i, 0) * rShapeDerivatives(j, 1));
            rOutput(FirstRow + 1, FirstCol + 1) +=
                DynamicViscosity *
                (FourThirds * rShapeDerivatives(i, 1) * rShapeDerivatives(j, 1) +
                 rShapeDerivatives(i, 0) * rShapeDerivatives(j, 0));

            // Update Counter
            FirstRow += 2;
        }
        FirstRow = 0;
        FirstCol += 2;
    }

    KRATOS_CATCH("");
}

template <>
void AddViscousTerm<3, 4>(BoundedMatrix<double, 12, 12>& rOutput,
                          const Matrix& rShapeDerivatives,
                          const double DynamicViscosity)
{
    KRATOS_TRY

    const IndexType number_of_nodes = 4;

    KRATOS_DEBUG_ERROR_IF(number_of_nodes != rShapeDerivatives.size1())
        << "Shape derivatives size1 and number of nodes mismatch [ "
        << rShapeDerivatives.size1() << " != " << number_of_nodes << " ].\n";

    KRATOS_DEBUG_ERROR_IF(rShapeDerivatives.size2() != 3)
        << "Shape derivatives size2 and dimension mismatch [ "
        << rShapeDerivatives.size2() << " != 3 ].\n";

    const double OneThird = 1.0 / 3.0;
    const double nTwoThirds = -2.0 / 3.0;

    IndexType FirstRow(0), FirstCol(0);

    for (IndexType j = 0; j < number_of_nodes; ++j)
    {
        for (IndexType i = 0; i < number_of_nodes; ++i)
        {
            // (dN_i/dx_k dN_j/dx_k)
            const double Diag = rShapeDerivatives(i, 0) * rShapeDerivatives(j, 0) +
                                rShapeDerivatives(i, 1) * rShapeDerivatives(j, 1) +
                                rShapeDerivatives(i, 2) * rShapeDerivatives(j, 2);

            // First Row
            rOutput(FirstRow, FirstCol) +=
                DynamicViscosity *
                (OneThird * rShapeDerivatives(i, 0) * rShapeDerivatives(j, 0) + Diag);
            rOutput(FirstRow, FirstCol + 1) +=
                DynamicViscosity *
                (nTwoThirds * rShapeDerivatives(i, 0) * rShapeDerivatives(j, 1) +
                 rShapeDerivatives(i, 1) * rShapeDerivatives(j, 0));
            rOutput(FirstRow, FirstCol + 2) +=
                DynamicViscosity *
                (nTwoThirds * rShapeDerivatives(i, 0) * rShapeDerivatives(j, 2) +
                 rShapeDerivatives(i, 2) * rShapeDerivatives(j, 0));

            // Second Row
            rOutput(FirstRow + 1, FirstCol) +=
                DynamicViscosity *
                (nTwoThirds * rShapeDerivatives(i, 1) * rShapeDerivatives(j, 0) +
                 rShapeDerivatives(i, 0) * rShapeDerivatives(j, 1));
            rOutput(FirstRow + 1, FirstCol + 1) +=
                DynamicViscosity *
                (OneThird * rShapeDerivatives(i, 1) * rShapeDerivatives(j, 1) + Diag);
            rOutput(FirstRow + 1, FirstCol + 2) +=
                DynamicViscosity *
                (nTwoThirds * rShapeDerivatives(i, 1) * rShapeDerivatives(j, 2) +
                 rShapeDerivatives(i, 2) * rShapeDerivatives(j, 1));

            // Third Row
            rOutput(FirstRow + 2, FirstCol) +=
                DynamicViscosity *
                (nTwoThirds * rShapeDerivatives(i, 2) * rShapeDerivatives(j, 0) +
                 rShapeDerivatives(i, 0) * rShapeDerivatives(j, 2));
            rOutput(FirstRow + 2, FirstCol + 1) +=
                DynamicViscosity *
                (nTwoThirds * rShapeDerivatives(i, 2) * rShapeDerivatives(j, 1) +
                 rShapeDerivatives(i, 1) * rShapeDerivatives(j, 2));
            rOutput(FirstRow + 2, FirstCol + 2) +=
                DynamicViscosity *
                (OneThird * rShapeDerivatives(i, 2) * rShapeDerivatives(j, 2) + Diag);

            // Update Counter
            FirstRow += 3;
        }
        FirstRow = 0;
        FirstCol += 3;
    }

    KRATOS_CATCH("");
}

template <IndexType TDim, IndexType TNumNodes>
BoundedMatrix<double, TDim * TNumNodes, TDim * TNumNodes> CalculateVelocityMatrixGaussPointContributions(
    const double Density,
    const double DynamicViscosity,
    const double TauOne,
    const double TauTwo,
    const BoundedVector<double, TNumNodes>& rConvectiveVelocity,
    const Vector& rShapeFunctions,
    const Matrix& rShapeDerivatives,
    const double GaussWeight)
{
    KRATOS_TRY

    constexpr IndexType local_size = TDim * TNumNodes;

    KRATOS_DEBUG_ERROR_IF(rShapeFunctions.size() != TNumNodes)
        << "Shape functions size and nodes size mismatch [ "
        << rShapeFunctions.size() << " != " << TNumNodes << " ].\n";

    KRATOS_DEBUG_ERROR_IF(rShapeDerivatives.size1() != TNumNodes)
        << "Shape function derviatives size1 and nodes size mismatch [ "
        << rShapeDerivatives.size1() << " != " << TNumNodes << " ].\n";

    KRATOS_DEBUG_ERROR_IF(rShapeDerivatives.size2() != TDim)
        << "Shape function derviatives size2 and dimension mismatch [ "
        << rShapeDerivatives.size2() << " != " << TDim << " ].\n";

    BoundedMatrix<double, local_size, local_size> output =
        ZeroMatrix(local_size, local_size);

    const double coeff_1 = Density * GaussWeight;
    const double coeff_2 = std::pow(Density, 2) * TauOne * GaussWeight;
    const double coeff_3 = TauTwo * GaussWeight;

    IndexType row_index{0}, col_index{0};
    double value{0.0};
    for (IndexType a = 0; a < TNumNodes; ++a)
    {
        row_index = a * TDim;
        for (IndexType b = 0; b < TNumNodes; ++b)
        {
            col_index = b * TDim;
            for (IndexType i = 0; i < TDim; ++i)
            {
                value = coeff_1 * rShapeFunctions[a] * rConvectiveVelocity[b];
                value += coeff_2 * rConvectiveVelocity[a] * rConvectiveVelocity[b];

                output(row_index + i, col_index + i) += value;

                for (IndexType j = 0; j < TDim; ++j)
                {
                    output(row_index + i, col_index + j) +=
                        coeff_3 * rShapeDerivatives(a, i) * rShapeDerivatives(b, j);
                }
            }
        }
    }

    AddViscousTerm<TDim, TNumNodes>(output, rShapeDerivatives, DynamicViscosity * GaussWeight);

    return output;

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
BoundedMatrix<double, TDim * TNumNodes, TNumNodes> CalculateVelocityPressureMatrixGaussPointContributions(
    const double Density,
    const double TauOne,
    const BoundedVector<double, TNumNodes>& rConvectiveVelocity,
    const Vector& rShapeFunctions,
    const Matrix& rShapeDerivatives,
    const double GaussWeight)
{
    KRATOS_TRY

    KRATOS_DEBUG_ERROR_IF(rShapeFunctions.size() != TNumNodes)
        << "Shape functions size and nodes size mismatch [ "
        << rShapeFunctions.size() << " != " << TNumNodes << " ].\n";

    KRATOS_DEBUG_ERROR_IF(rShapeDerivatives.size1() != TNumNodes)
        << "Shape function derviatives size1 and nodes size mismatch [ "
        << rShapeDerivatives.size1() << " != " << TNumNodes << " ].\n";

    KRATOS_DEBUG_ERROR_IF(rShapeDerivatives.size2() != TDim)
        << "Shape function derviatives size2 and dimension mismatch [ "
        << rShapeDerivatives.size2() << " != " << TDim << " ].\n";

    BoundedMatrix<double, TDim * TNumNodes, TNumNodes> output;

    const double coeff_1 = -1.0 * GaussWeight;
    const double coeff_2 = Density * TauOne * GaussWeight;

    IndexType row_index{0};
    double value{0.0};
    for (IndexType a = 0; a < TNumNodes; ++a)
    {
        row_index = a * TDim;
        for (IndexType b = 0; b < TNumNodes; ++b)
        {
            for (IndexType i = 0; i < TDim; ++i)
            {
                value = coeff_1 * rShapeDerivatives(a, i) * rShapeFunctions[b];
                value += coeff_2 * rConvectiveVelocity[a] * rShapeDerivatives(b, i);

                output(row_index + i, b) = value;
            }
        }
    }

    return output;

    KRATOS_CATCH("");
}

template <>
BoundedVector<double, 2> GetValuesVector<double, 2>(const GeometryType& rGeometry,
                                                    const Variable<double>& rVariable,
                                                    const int Step)
{
    BoundedVector<double, 2> output;
    output[0] = rGeometry[0].FastGetSolutionStepValue(rVariable, Step);
    output[1] = rGeometry[1].FastGetSolutionStepValue(rVariable, Step);
    return output;
}

template <>
BoundedVector<double, 3> GetValuesVector<double, 3>(const GeometryType& rGeometry,
                                                    const Variable<double>& rVariable,
                                                    const int Step)
{
    BoundedVector<double, 3> output;
    output[0] = rGeometry[0].FastGetSolutionStepValue(rVariable, Step);
    output[1] = rGeometry[1].FastGetSolutionStepValue(rVariable, Step);
    output[2] = rGeometry[2].FastGetSolutionStepValue(rVariable, Step);
    return output;
}

template <>
BoundedVector<double, 4> GetValuesVector<double, 4>(const GeometryType& rGeometry,
                                                    const Variable<double>& rVariable,
                                                    const int Step)
{
    BoundedVector<double, 4> output;
    output[0] = rGeometry[0].FastGetSolutionStepValue(rVariable, Step);
    output[1] = rGeometry[1].FastGetSolutionStepValue(rVariable, Step);
    output[2] = rGeometry[2].FastGetSolutionStepValue(rVariable, Step);
    output[3] = rGeometry[3].FastGetSolutionStepValue(rVariable, Step);
    return output;
}

template <>
BoundedVector<double, 4> GetValuesVector<array_1d<double, 3>, 4>(
    const GeometryType& rGeometry, const Variable<array_1d<double, 3>>& rVariable, const int Step)
{
    BoundedVector<double, 4> output;

    IndexType local_index{0};

    for (IndexType a = 0; a < 2; ++a)
    {
        const array_1d<double, 3>& r_value =
            rGeometry[a].FastGetSolutionStepValue(rVariable, Step);
        output[local_index++] = r_value[0];
        output[local_index++] = r_value[1];
    }

    return output;
}

template <>
BoundedVector<double, 6> GetValuesVector<array_1d<double, 3>, 6>(
    const GeometryType& rGeometry, const Variable<array_1d<double, 3>>& rVariable, const int Step)
{
    BoundedVector<double, 6> output;

    IndexType local_index{0};

    for (IndexType a = 0; a < 3; ++a)
    {
        const array_1d<double, 3>& r_value =
            rGeometry[a].FastGetSolutionStepValue(rVariable, Step);
        output[local_index++] = r_value[0];
        output[local_index++] = r_value[1];
    }

    return output;
}

template <>
BoundedVector<double, 9> GetValuesVector<array_1d<double, 3>, 9>(
    const GeometryType& rGeometry, const Variable<array_1d<double, 3>>& rVariable, const int Step)
{
    BoundedVector<double, 9> output;

    IndexType local_index{0};

    for (IndexType a = 0; a < 3; ++a)
    {
        const array_1d<double, 3>& r_value =
            rGeometry[a].FastGetSolutionStepValue(rVariable, Step);
        output[local_index++] = r_value[0];
        output[local_index++] = r_value[1];
        output[local_index++] = r_value[2];
    }

    return output;
}

template <>
BoundedVector<double, 12> GetValuesVector<array_1d<double, 3>, 12>(
    const GeometryType& rGeometry, const Variable<array_1d<double, 3>>& rVariable, const int Step)
{
    BoundedVector<double, 12> output;

    IndexType local_index{0};

    for (IndexType a = 0; a < 4; ++a)
    {
        const array_1d<double, 3>& r_value =
            rGeometry[a].FastGetSolutionStepValue(rVariable, Step);
        output[local_index++] = r_value[0];
        output[local_index++] = r_value[1];
        output[local_index++] = r_value[2];
    }

    return output;
}

// template instantiations

template BoundedVector<double, 3> GetConvectionOperator<2, 3>(const array_1d<double, 3>&,
                                                              const Matrix&);

template BoundedVector<double, 4> GetConvectionOperator<3, 4>(const array_1d<double, 3>&,
                                                              const Matrix&);

template BoundedMatrix<double, 6, 6> CalculateVelocityMatrixGaussPointContributions<2, 3>(
    const double,
    const double,
    const double,
    const double,
    const BoundedVector<double, 3>&,
    const Vector&,
    const Matrix&,
    const double);

template BoundedMatrix<double, 12, 12> CalculateVelocityMatrixGaussPointContributions<3, 4>(
    const double,
    const double,
    const double,
    const double,
    const BoundedVector<double, 4>&,
    const Vector&,
    const Matrix&,
    const double);

template BoundedMatrix<double, 6, 3> CalculateVelocityPressureMatrixGaussPointContributions<2, 3>(
    const double, const double, const BoundedVector<double, 3>&, const Vector&, const Matrix&, const double);

template BoundedMatrix<double, 12, 4> CalculateVelocityPressureMatrixGaussPointContributions<3, 4>(
    const double, const double, const BoundedVector<double, 4>&, const Vector&, const Matrix&, const double);

} // namespace SegregatedVMSElementUtilities
///@}

} // namespace Kratos