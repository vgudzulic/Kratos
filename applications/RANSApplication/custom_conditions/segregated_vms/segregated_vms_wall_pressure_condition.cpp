//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <string>

// External includes

// Project includes
#include "includes/cfd_variables.h"
#include "includes/checks.h"
#include "includes/define.h"
#include "includes/variables.h"

// Application includes
#include "custom_utilities/rans_calculation_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "segregated_vms_wall_pressure_condition.h"

namespace Kratos
{
template <unsigned int TDim, unsigned int TNumNodes>
void SegregatedVMSWallPressureCondition<TDim, TNumNodes>::EquationIdVector(
    EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
{
    if (rResult.size() != TNumNodes)
        rResult.resize(TNumNodes, false);

    for (unsigned int i_node = 0; i_node < TNumNodes; ++i_node)
    {
        rResult[i_node] = this->GetGeometry()[i_node].GetDof(PRESSURE).EquationId();
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void SegregatedVMSWallPressureCondition<TDim, TNumNodes>::GetDofList(DofsVectorType& rElementalDofList,
                                                                     ProcessInfo& rCurrentProcessInfo)
{
    if (rElementalDofList.size() != TNumNodes)
        rElementalDofList.resize(TNumNodes);

    for (unsigned int i_node = 0; i_node < TNumNodes; ++i_node)
    {
        rElementalDofList[i_node] = this->GetGeometry()[i_node].pGetDof(PRESSURE);
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void SegregatedVMSWallPressureCondition<TDim, TNumNodes>::GetFirstDerivativesVector(Vector& Values, int Step)
{
    if (Values.size() != TNumNodes)
        Values.resize(TNumNodes, false);

    for (unsigned int i_node = 0; i_node < 2; ++i_node)
    {
        Values[i_node] =
            this->GetGeometry()[i_node].FastGetSolutionStepValue(PRESSURE, Step);
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void SegregatedVMSWallPressureCondition<TDim, TNumNodes>::GetSecondDerivativesVector(Vector& Values, int Step)
{
    if (Values.size() != TNumNodes)
        Values.resize(TNumNodes, false);

    noalias(Values) = ZeroVector(TNumNodes);
}

template <unsigned int TDim, unsigned int TNumNodes>
int SegregatedVMSWallPressureCondition<TDim, TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    int check = BaseType::Check(rCurrentProcessInfo);

    const GeometryType& r_geometry = this->GetGeometry();

    for (IndexType i_node = 0; i_node < TNumNodes; ++i_node)
    {
        const NodeType& r_node = r_geometry[i_node];

        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(MESH_VELOCITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DENSITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(BODY_FORCE, r_node);

        KRATOS_CHECK_DOF_IN_NODE(PRESSURE, r_node);
    }

    return check;

    KRATOS_CATCH("");
}

template <>
void SegregatedVMSWallPressureCondition<2, 2>::CalculateNormal()
{
    KRATOS_TRY

    Geometry<Node<3>>& pGeometry = this->GetGeometry();

    array_1d<double, 3> normal;

    normal[0] = pGeometry[1].Y() - pGeometry[0].Y();
    normal[1] = -(pGeometry[1].X() - pGeometry[0].X());
    normal[2] = 0.00;

    const double normal_magnitude = norm_2(normal);
    KRATOS_ERROR_IF(normal_magnitude == 0.0)
        << "Normal magnitude is zero in element " << this->Info();

    noalias(normal) = normal * (1 / normal_magnitude);
    this->SetValue(NORMAL, normal);

    KRATOS_CATCH("");
}

template <>
void SegregatedVMSWallPressureCondition<3, 3>::CalculateNormal()
{
    Geometry<Node<3>>& pGeometry = this->GetGeometry();

    array_1d<double, 3> v1, v2, normal;

    v1[0] = pGeometry[1].X() - pGeometry[0].X();
    v1[1] = pGeometry[1].Y() - pGeometry[0].Y();
    v1[2] = pGeometry[1].Z() - pGeometry[0].Z();

    v2[0] = pGeometry[2].X() - pGeometry[0].X();
    v2[1] = pGeometry[2].Y() - pGeometry[0].Y();
    v2[2] = pGeometry[2].Z() - pGeometry[0].Z();

    MathUtils<double>::CrossProduct(normal, v1, v2);

    const double normal_magnitude = norm_2(normal);
    KRATOS_ERROR_IF(normal_magnitude == 0.0)
        << "Normal magnitude is zero in element " << this->Info();

    noalias(normal) = normal * (1 / normal_magnitude);
    this->SetValue(NORMAL, normal);
}

template <unsigned int TDim, unsigned int TNumNodes>
void SegregatedVMSWallPressureCondition<TDim, TNumNodes>::Initialize()
{
    KRATOS_TRY;

    if (RansCalculationUtilities::IsInlet(*this))
    {
        const ConditionType& r_parent_condition =
            *(this->GetValue(PARENT_CONDITION_POINTER));

        this->CalculateNormal();

        KRATOS_ERROR_IF(r_parent_condition.GetValue(NEIGHBOUR_ELEMENTS).size() == 0)
            << this->Info() << " cannot find parent element\n";
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
std::string SegregatedVMSWallPressureCondition<TDim, TNumNodes>::Info() const
{
    std::stringstream buffer;
    buffer << "SegregatedVMSWallPressureCondition" << TDim << "D";
    return buffer.str();
}

template <unsigned int TDim, unsigned int TNumNodes>
void SegregatedVMSWallPressureCondition<TDim, TNumNodes>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "SegregatedVMSWallPressureCondition";
}

template <unsigned int TDim, unsigned int TNumNodes>
void SegregatedVMSWallPressureCondition<TDim, TNumNodes>::PrintData(std::ostream& rOStream) const
{
}

template <unsigned int TDim, unsigned int TNumNodes>
void SegregatedVMSWallPressureCondition<TDim, TNumNodes>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    if (rLeftHandSideMatrix.size1() != TNumNodes)
        rLeftHandSideMatrix.resize(TNumNodes, TNumNodes);

    if (rRightHandSideVector.size() != TNumNodes)
        rRightHandSideVector.resize(TNumNodes);

    noalias(rLeftHandSideMatrix) = ZeroMatrix(TNumNodes, TNumNodes);
    noalias(rRightHandSideVector) = ZeroVector(TNumNodes);
}

template <unsigned int TDim, unsigned int TNumNodes>
void SegregatedVMSWallPressureCondition<TDim, TNumNodes>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{
    if (rLeftHandSideMatrix.size1() != TNumNodes)
        rLeftHandSideMatrix.resize(TNumNodes, TNumNodes);

    noalias(rLeftHandSideMatrix) = ZeroMatrix(TNumNodes, TNumNodes);
}

template <unsigned int TDim, unsigned int TNumNodes>
void SegregatedVMSWallPressureCondition<TDim, TNumNodes>::CalculateRightHandSide(
    VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    if (rRightHandSideVector.size() != TNumNodes)
        rRightHandSideVector.resize(TNumNodes);

    noalias(rRightHandSideVector) = ZeroVector(TNumNodes);

    // if (RansCalculationUtilities::IsInlet(*this))
    // {
    //     // calculates parent element pressure gradient
    //     const ConditionType& r_parent_condition =
    //         *(this->GetValue(PARENT_CONDITION_POINTER));
    //     const ModelPart::ElementType& r_parent_element =
    //         r_parent_condition.GetValue(NEIGHBOUR_ELEMENTS)[0];
    //     const GeometryType& r_parent_geometry = r_parent_element.GetGeometry();

    //     // Get Shape function data
    //     Vector parent_element_gauss_weights;
    //     Matrix parent_element_shape_functions;
    //     ShapeFunctionDerivativesArrayType parent_element_shape_derivatives;
    //     RansCalculationUtilities::CalculateGeometryData(
    //         r_parent_geometry, GeometryData::GI_GAUSS_1, parent_element_gauss_weights,
    //         parent_element_shape_functions, parent_element_shape_derivatives);

    //     array_1d<double, 3> pressure_gradient;
    //     const Matrix& r_parent_element_shape_derivatives =
    //         parent_element_shape_derivatives[0];
    //     const Vector parent_elemnet_gauss_shape_functions =
    //         row(parent_element_shape_functions, 0);

    //     const double density = RansCalculationUtilities::EvaluateInPoint(
    //         r_parent_geometry, PRESSURE, parent_elemnet_gauss_shape_functions);
    //     const array_1d<double, 3>& r_body_force =
    //         RansCalculationUtilities::EvaluateInPoint(
    //             r_parent_geometry, BODY_FORCE, parent_elemnet_gauss_shape_functions) *
    //         density;
    //     const array_1d<double, 3>& r_velocity = RansCalculationUtilities::EvaluateInPoint(
    //         r_parent_geometry, VELOCITY, parent_elemnet_gauss_shape_functions);
    //     BoundedMatrix<double, TDim, TDim> velocity_gradient;
    //     RansCalculationUtilities::CalculateGradient<TDim>(
    //         velocity_gradient, r_parent_geometry, VELOCITY, r_parent_element_shape_derivatives);

    //     noalias(pressure_gradient) = r_body_force;
    //     for (IndexType i = 0; i < TDim; ++i)
    //     {
    //         for (IndexType j = 0; j < TDim; ++j)
    //         {
    //             pressure_gradient[i] -=
    //                 density * r_velocity[j] * velocity_gradient(i, j);
    //         }
    //     }

    //     const array_1d<double, 3>& r_normal = this->GetValue(NORMAL);
    //     const double pressure_gradient_normal = inner_prod(pressure_gradient, r_normal);

    //     // Get Shape function data
    //     MatrixType shape_functions;
    //     VectorType gauss_weights;
    //     RansCalculationUtilities::CalculateConditionGeometryData(
    //         this->GetGeometry(), GeometryData::GI_GAUSS_2, gauss_weights, shape_functions);
    //     const IndexType number_of_gauss_points = gauss_weights.size();

    //     for (IndexType g = 0; g < number_of_gauss_points; ++g)
    //     {
    //         const Vector& gauss_shape_functions = row(shape_functions, g);
    //         const double coeff_1 = gauss_weights[g] * pressure_gradient_normal;

    //         for (IndexType a = 0; a < TNumNodes; ++a)
    //         {
    //             rRightHandSideVector[a] += gauss_shape_functions[a] * coeff_1;
    //         }
    //     }
    // }
}

template <unsigned int TDim, unsigned int TNumNodes>
void SegregatedVMSWallPressureCondition<TDim, TNumNodes>::CalculateDampingMatrix(
    MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo)
{
    VectorType RHS;
    this->CalculateLocalVelocityContribution(rDampingMatrix, RHS, rCurrentProcessInfo);
}

template <unsigned int TDim, unsigned int TNumNodes>
void SegregatedVMSWallPressureCondition<TDim, TNumNodes>::CalculateLocalVelocityContribution(
    MatrixType& rDampMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    if (rDampMatrix.size1() != TNumNodes)
        rDampMatrix.resize(TNumNodes, TNumNodes);
    if (rRightHandSideVector.size() != TNumNodes)
        rRightHandSideVector.resize(TNumNodes);

    noalias(rDampMatrix) = ZeroMatrix(TNumNodes, TNumNodes);
}

template <unsigned int TDim, unsigned int TNumNodes>
void SegregatedVMSWallPressureCondition<TDim, TNumNodes>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
}

template <unsigned int TDim, unsigned int TNumNodes>
void SegregatedVMSWallPressureCondition<TDim, TNumNodes>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
}

// template instantiations

template class SegregatedVMSWallPressureCondition<2, 2>;
template class SegregatedVMSWallPressureCondition<3, 3>;

} // namespace Kratos.
