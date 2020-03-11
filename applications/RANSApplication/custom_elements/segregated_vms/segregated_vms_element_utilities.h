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

#if !defined(KRATOS_SEGREGATED_VMS_ELEMENT_UTILITIES_H_INCLUDED)
#define KRATOS_SEGREGATED_VMS_ELEMENT_UTILITIES_H_INCLUDED

// System includes
#include <cmath>

// Project includes
#include "geometries/geometry.h"
#include "geometries/geometry_data.h"
#include "includes/model_part.h"

namespace Kratos
{
///@name Kratos Globals
///@{

namespace SegregatedVMSElementUtilities
{
/// Node type
using NodeType = ModelPart::NodeType;
using ElementType = ModelPart::ElementType;
using ConditionType = ModelPart::ConditionType;

/// Geometry type (using with given NodeType)
using GeometryType = Geometry<NodeType>;

using IndexType = unsigned int;
using SizeType = std::size_t;

template <unsigned int TDim, unsigned int TNumNodes>
BoundedVector<double, TNumNodes> GetConvectionOperator(const array_1d<double, 3>& rVector,
                                                       const Matrix& rShapeDerivatives);

double CalculateTauOne(const double VelocityMagnitude,
                       const double ElemSize,
                       const double Density,
                       const double Viscosity,
                       const double DynamicTau,
                       const double DeltaTime);

double CalculateTauTwo(const double VelocityMagnitude,
                       const double ElemSize,
                       const double Density,
                       const double Viscosity);

template <unsigned int TDim, unsigned int TNumNodes>
void AddViscousTerm(BoundedMatrix<double, TDim * TNumNodes, TDim * TNumNodes>& rOutput,
                    const Matrix& rShapeDerivatives,
                    const double DynamicViscosity);

template <unsigned int TDim, unsigned int TNumNodes>
BoundedMatrix<double, TDim * TNumNodes, TDim * TNumNodes> CalculateVelocityMatrixGaussPointContributions(
    const double Density,
    const double DynamicViscosity,
    const double TauOne,
    const double TauTwo,
    const BoundedVector<double, TNumNodes>& rConvectiveVelocity,
    const Vector& rShapeFunctions,
    const Matrix& rShapeDerivatives,
    const double GaussWeight);

template <unsigned int TDim, unsigned int TNumNodes>
BoundedMatrix<double, TDim * TNumNodes, TNumNodes> CalculateVelocityPressureMatrixGaussPointContributions(
    const double Density,
    const double TauOne,
    const BoundedVector<double, TNumNodes>& rConvectiveVelocity,
    const Vector& rShapeFunctions,
    const Matrix& rShapeDerivatives,
    const double GaussWeight);

template <typename TDataType, unsigned int TSize>
BoundedVector<double, TSize> GetValuesVector(const GeometryType& rGeometry,
                                             const Variable<TDataType>& rVariable,
                                             const int Step = 0);

} // namespace SegregatedVMSElementUtilities

///@}

} // namespace Kratos

#endif // KRATOS_SEGREGATED_VMS_ELEMENT_UTILITIES_H_INCLUDED defined