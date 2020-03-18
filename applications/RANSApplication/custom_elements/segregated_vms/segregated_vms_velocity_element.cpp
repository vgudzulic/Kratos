//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes

// External includes

// Project includes
#include "includes/cfd_variables.h"
#include "includes/checks.h"
#include "includes/variables.h"

// Application includes
#include "custom_utilities/rans_calculation_utilities.h"
#include "segregated_vms_element_utilities.h"

// Include base h
#include "segregated_vms_velocity_element.h"

namespace Kratos
{
template <unsigned int TDim, unsigned int TNumNodes>
GeometryData::IntegrationMethod SegregatedVMSVelocityElement<TDim, TNumNodes>::GetIntegrationMethod() const
{
    return GeometryData::GI_GAUSS_1;
}

template <unsigned int TDim, unsigned int TNumNodes>
void SegregatedVMSVelocityElement<TDim, TNumNodes>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    // Check sizes and initialize matrix
    if (rLeftHandSideMatrix.size1() != LocalSize || rLeftHandSideMatrix.size2() != LocalSize)
        rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);

    noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize);

    // Calculate RHS
    this->CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
}

template <unsigned int TDim, unsigned int TNumNodes>
void SegregatedVMSVelocityElement<TDim, TNumNodes>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{
    // Check sizes and initialize matrix
    if (rLeftHandSideMatrix.size1() != LocalSize || rLeftHandSideMatrix.size2() != LocalSize)
        rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);

    noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize);
}

template <unsigned int TDim, unsigned int TNumNodes>
void SegregatedVMSVelocityElement<TDim, TNumNodes>::CalculateMassMatrix(MatrixType& rMassMatrix,
                                                                        ProcessInfo& rCurrentProcessInfo)
{
    // Check sizes and initialize matrix
    if (rMassMatrix.size1() != LocalSize || rMassMatrix.size2() != LocalSize)
        rMassMatrix.resize(LocalSize, LocalSize, false);

    noalias(rMassMatrix) = ZeroMatrix(LocalSize, LocalSize);
}

template <unsigned int TDim, unsigned int TNumNodes>
void SegregatedVMSVelocityElement<TDim, TNumNodes>::CalculateGeometryData(
    Vector& rGaussWeights, Matrix& rNContainer, GeometryType::ShapeFunctionsGradientsType& rDN_DX) const
{
    RansCalculationUtilities::CalculateGeometryData(
        this->GetGeometry(), this->GetIntegrationMethod(), rGaussWeights, rNContainer, rDN_DX);
}

template <unsigned int TDim, unsigned int TNumNodes>
double SegregatedVMSVelocityElement<TDim, TNumNodes>::EvaluateInPoint(
    const Variable<double>& rVariable, const Vector& rShapeFunction, const int Step) const
{
    return RansCalculationUtilities::EvaluateInPoint(
        this->GetGeometry(), rVariable, rShapeFunction, Step);
}

template <unsigned int TDim, unsigned int TNumNodes>
array_1d<double, 3> SegregatedVMSVelocityElement<TDim, TNumNodes>::EvaluateInPoint(
    const Variable<array_1d<double, 3>>& rVariable, const Vector& rShapeFunction, const int Step) const
{
    return RansCalculationUtilities::EvaluateInPoint(
        this->GetGeometry(), rVariable, rShapeFunction, Step);
}

template <unsigned int TDim, unsigned int TNumNodes>
BoundedVector<double, TNumNodes> SegregatedVMSVelocityElement<TDim, TNumNodes>::GetConvectionOperator(
    const array_1d<double, 3>& rVector, const Matrix& rShapeDerivatives) const
{
    return SegregatedVMSElementUtilities::GetConvectionOperator<TDim, TNumNodes>(
        rVector, rShapeDerivatives);
}

template <unsigned int TDim, unsigned int TNumNodes>
array_1d<double, 3> SegregatedVMSVelocityElement<TDim, TNumNodes>::GetVelocity(
    const Vector& rShapeFunction, const int Step)
{
    const array_1d<double, 3>& r_velocity =
        this->EvaluateInPoint(VELOCITY, rShapeFunction, Step);
    const array_1d<double, 3>& r_mesh_velocity =
        this->EvaluateInPoint(MESH_VELOCITY, rShapeFunction, Step);
    return (r_velocity - r_mesh_velocity);
}

template <unsigned int TDim, unsigned int TNumNodes>
void SegregatedVMSVelocityElement<TDim, TNumNodes>::GetFirstDerivativesVector(Vector& Values, int Step)
{
    this->GetValuesVector(Values, Step);
}

template <>
void SegregatedVMSVelocityElement<2, 3>::EquationIdVector(EquationIdVectorType& rResult,
                                                          ProcessInfo& rCurrentProcessInfo)
{
    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize);

    IndexType local_index = 0;
    IndexType velocity_index = this->GetGeometry()[0].GetDofPosition(VELOCITY_X);

    for (IndexType i_node = 0; i_node < 3; ++i_node)
    {
        rResult[local_index++] =
            this->GetGeometry()[i_node].GetDof(VELOCITY_X, velocity_index).EquationId();
        rResult[local_index++] =
            this->GetGeometry()[i_node].GetDof(VELOCITY_Y, velocity_index + 1).EquationId();
    }
}

template <>
void SegregatedVMSVelocityElement<3, 4>::EquationIdVector(EquationIdVectorType& rResult,
                                                          ProcessInfo& rCurrentProcessInfo)
{
    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize);

    IndexType local_index = 0;
    IndexType velocity_index = this->GetGeometry()[0].GetDofPosition(VELOCITY_X);

    for (IndexType i_node = 0; i_node < 4; ++i_node)
    {
        rResult[local_index++] =
            this->GetGeometry()[i_node].GetDof(VELOCITY_X, velocity_index).EquationId();
        rResult[local_index++] =
            this->GetGeometry()[i_node].GetDof(VELOCITY_Y, velocity_index + 1).EquationId();
        rResult[local_index++] =
            this->GetGeometry()[i_node].GetDof(VELOCITY_Z, velocity_index + 2).EquationId();
    }
}

template <>
void SegregatedVMSVelocityElement<2, 3>::GetDofList(DofsVectorType& rElementalDofList,
                                                    ProcessInfo& rCurrentProcessInfo)
{
    if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    IndexType local_index = 0;

    for (IndexType i_node = 0; i_node < 3; ++i_node)
    {
        rElementalDofList[local_index++] = this->GetGeometry()[i_node].pGetDof(VELOCITY_X);
        rElementalDofList[local_index++] = this->GetGeometry()[i_node].pGetDof(VELOCITY_Y);
    }
}

template <>
void SegregatedVMSVelocityElement<3, 4>::GetDofList(DofsVectorType& rElementalDofList,
                                                    ProcessInfo& rCurrentProcessInfo)
{
    if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    IndexType local_index = 0;

    for (IndexType i_node = 0; i_node < 4; ++i_node)
    {
        rElementalDofList[local_index++] = this->GetGeometry()[i_node].pGetDof(VELOCITY_X);
        rElementalDofList[local_index++] = this->GetGeometry()[i_node].pGetDof(VELOCITY_Y);
        rElementalDofList[local_index++] = this->GetGeometry()[i_node].pGetDof(VELOCITY_Z);
    }
}

template <>
void SegregatedVMSVelocityElement<2, 3>::GetValuesVector(Vector& Values, int Step)
{
    if (Values.size() != LocalSize)
        Values.resize(LocalSize, false);

    IndexType local_index = 0;
    for (IndexType i_node = 0; i_node < 3; ++i_node)
    {
        const array_1d<double, 3>& r_velocity =
            this->GetGeometry()[i_node].FastGetSolutionStepValue(VELOCITY, Step);
        Values[local_index++] = r_velocity[0];
        Values[local_index++] = r_velocity[1];
    }
}

template <>
void SegregatedVMSVelocityElement<3, 4>::GetValuesVector(Vector& Values, int Step)
{
    if (Values.size() != LocalSize)
        Values.resize(LocalSize, false);

    IndexType local_index = 0;
    for (IndexType i_node = 0; i_node < 3; ++i_node)
    {
        const array_1d<double, 3>& r_velocity =
            this->GetGeometry()[i_node].FastGetSolutionStepValue(VELOCITY, Step);
        Values[local_index++] = r_velocity[0];
        Values[local_index++] = r_velocity[1];
        Values[local_index++] = r_velocity[2];
    }
}

template <>
void SegregatedVMSVelocityElement<2, 3>::GetSecondDerivativesVector(Vector& Values, int Step)
{
    if (Values.size() != LocalSize)
        Values.resize(LocalSize, false);

    IndexType local_index = 0;
    for (IndexType i_node = 0; i_node < 3; ++i_node)
    {
        const array_1d<double, 3>& r_acceleration =
            this->GetGeometry()[i_node].FastGetSolutionStepValue(ACCELERATION, Step);
        Values[local_index++] = r_acceleration[0];
        Values[local_index++] = r_acceleration[1];
    }
}

template <>
void SegregatedVMSVelocityElement<3, 4>::GetSecondDerivativesVector(Vector& Values, int Step)
{
    if (Values.size() != LocalSize)
        Values.resize(LocalSize, false);

    IndexType local_index = 0;
    for (IndexType i_node = 0; i_node < 4; ++i_node)
    {
        const array_1d<double, 3>& r_acceleration =
            this->GetGeometry()[i_node].FastGetSolutionStepValue(ACCELERATION, Step);
        Values[local_index++] = r_acceleration[0];
        Values[local_index++] = r_acceleration[1];
        Values[local_index++] = r_acceleration[2];
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
int SegregatedVMSVelocityElement<TDim, TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    int value = BaseType::Check(rCurrentProcessInfo);

    for (IndexType i = 0; i < this->GetGeometry().size(); ++i)
    {
        const NodeType& r_node = this->GetGeometry()[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(MESH_VELOCITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ACCELERATION, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DENSITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VISCOSITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(BODY_FORCE, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(PRESSURE, r_node);

        KRATOS_CHECK_DOF_IN_NODE(VELOCITY_X, r_node);
        KRATOS_CHECK_DOF_IN_NODE(VELOCITY_Y, r_node);
        if (TDim == 3)
            KRATOS_CHECK_DOF_IN_NODE(VELOCITY_Z, r_node);
    }

    return value;

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void SegregatedVMSVelocityElement<TDim, TNumNodes>::CalculateRightHandSide(
    VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Check sizes and initialize
    if (rRightHandSideVector.size() != LocalSize)
        rRightHandSideVector.resize(LocalSize, false);

    noalias(rRightHandSideVector) = ZeroVector(LocalSize);

    // Get Shape function data
    Vector gauss_weights;
    Matrix shape_functions;
    ShapeFunctionDerivativesArrayType shape_derivatives;
    this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
    const IndexType num_gauss_points = gauss_weights.size();

    double value = 0.0;
    for (IndexType g = 0; g < num_gauss_points; ++g)
    {
        const Matrix& r_shape_derivatives = shape_derivatives[g];
        const Vector& gauss_shape_functions = row(shape_functions, g);
        const double weight = gauss_weights[g];

        const double density = this->EvaluateInPoint(DENSITY, gauss_shape_functions);
        const double pressure = this->EvaluateInPoint(PRESSURE, gauss_shape_functions);
        const array_1d<double, 3>& r_body_force =
            this->EvaluateInPoint(BODY_FORCE, gauss_shape_functions) * density;

        for (IndexType a = 0; a < TNumNodes; ++a)
        {
            IndexType row_index = a * TDim;
            for (IndexType d = 0; d < TDim; ++d)
            {
                value = gauss_shape_functions[a] * r_body_force[d] * weight;
                value += r_shape_derivatives(a, d) * pressure * weight;
                rRightHandSideVector[row_index + d] += value;
            }
        }
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void SegregatedVMSVelocityElement<TDim, TNumNodes>::CalculateLocalVelocityContribution(
    MatrixType& rDampingMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Check sizes and initialize matrix
    if (rDampingMatrix.size1() != LocalSize || rDampingMatrix.size2() != LocalSize)
        rDampingMatrix.resize(LocalSize, LocalSize, false);

    noalias(rDampingMatrix) = ZeroMatrix(LocalSize, LocalSize);

    // Get Shape function data
    Vector gauss_weights;
    Matrix shape_functions;
    ShapeFunctionDerivativesArrayType shape_derivatives;
    this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
    const IndexType num_gauss_points = gauss_weights.size();

    const double dynamic_tau = rCurrentProcessInfo[DYNAMIC_TAU];
    const double delta_time = rCurrentProcessInfo[DELTA_TIME];
    const double element_length = this->GetGeometry().Length();

    for (IndexType g = 0; g < num_gauss_points; ++g)
    {
        const Matrix& r_shape_derivatives = shape_derivatives[g];
        const Vector& gauss_shape_functions = row(shape_functions, g);

        const double density = this->EvaluateInPoint(DENSITY, gauss_shape_functions);
        const double dynamic_viscosity =
            this->EvaluateInPoint(VISCOSITY, gauss_shape_functions) * density;

        const array_1d<double, 3>& r_velocity = this->GetVelocity(gauss_shape_functions);
        const double velocity_magnitude = norm_2(r_velocity);
        const BoundedVector<double, TNumNodes>& r_convective_velocity =
            this->GetConvectionOperator(r_velocity, r_shape_derivatives);

        const double tau_one = SegregatedVMSElementUtilities::CalculateTauOne(
            velocity_magnitude, element_length, density, dynamic_viscosity,
            dynamic_tau, delta_time);
        const double tau_two = SegregatedVMSElementUtilities::CalculateTauTwo(
            velocity_magnitude, element_length, density, dynamic_viscosity);

        noalias(rDampingMatrix) +=
            SegregatedVMSElementUtilities::CalculateVelocityMatrixGaussPointContributions<TDim, TNumNodes>(
                density, dynamic_viscosity, tau_one, tau_two, r_convective_velocity,
                gauss_shape_functions, r_shape_derivatives, gauss_weights[g]);

        const array_1d<double, 3>& r_body_force =
            this->EvaluateInPoint(BODY_FORCE, gauss_shape_functions) * density;

        array_1d<double, 3> pressure_gradient;
        RansCalculationUtilities::CalculateGradient(
            pressure_gradient, this->GetGeometry(), PRESSURE, r_shape_derivatives);

        noalias(pressure_gradient) -= r_body_force;

        const double coeff_1 = density * tau_one * gauss_weights[g];

        IndexType row_index{0};
        for (IndexType a = 0; a < TNumNodes; ++a)
        {
            row_index = a * TDim;
            for (IndexType b = 0; b < TNumNodes; ++b)
            {
                for (IndexType i = 0; i < TDim; ++i)
                {
                    rRightHandSideVector[row_index + i] -=
                        coeff_1 * r_convective_velocity[a] * pressure_gradient[i];
                }
            }
        }
    }

    VectorType values;
    this->GetValuesVector(values);
    noalias(rRightHandSideVector) -= prod(rDampingMatrix, values);

    KRATOS_CATCH("");
}

// template instantiations

template class SegregatedVMSVelocityElement<2>;
template class SegregatedVMSVelocityElement<3>;

} // namespace Kratos.
