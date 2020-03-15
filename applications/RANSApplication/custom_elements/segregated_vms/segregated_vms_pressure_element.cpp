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
#include "segregated_vms_pressure_element.h"

namespace Kratos
{
template <unsigned int TDim, unsigned int TNumNodes>
GeometryData::IntegrationMethod SegregatedVMSPressureElement<TDim, TNumNodes>::GetIntegrationMethod() const
{
    return GeometryData::GI_GAUSS_1;
}

template <unsigned int TDim, unsigned int TNumNodes>
void SegregatedVMSPressureElement<TDim, TNumNodes>::CalculateLocalSystem(
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
void SegregatedVMSPressureElement<TDim, TNumNodes>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{
    // Check sizes and initialize matrix
    if (rLeftHandSideMatrix.size1() != LocalSize || rLeftHandSideMatrix.size2() != LocalSize)
        rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);

    noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize);
}

template <unsigned int TDim, unsigned int TNumNodes>
void SegregatedVMSPressureElement<TDim, TNumNodes>::CalculateMassMatrix(MatrixType& rMassMatrix,
                                                                        ProcessInfo& rCurrentProcessInfo)
{
    // Check sizes and initialize matrix
    if (rMassMatrix.size1() != LocalSize || rMassMatrix.size2() != LocalSize)
        rMassMatrix.resize(LocalSize, LocalSize, false);

    noalias(rMassMatrix) = ZeroMatrix(LocalSize, LocalSize);
}

template <unsigned int TDim, unsigned int TNumNodes>
void SegregatedVMSPressureElement<TDim, TNumNodes>::CalculateGeometryData(
    Vector& rGaussWeights, Matrix& rNContainer, GeometryType::ShapeFunctionsGradientsType& rDN_DX) const
{
    RansCalculationUtilities::CalculateGeometryData(
        this->GetGeometry(), this->GetIntegrationMethod(), rGaussWeights, rNContainer, rDN_DX);
}

template <unsigned int TDim, unsigned int TNumNodes>
double SegregatedVMSPressureElement<TDim, TNumNodes>::EvaluateInPoint(
    const Variable<double>& rVariable, const Vector& rShapeFunction, const int Step) const
{
    return RansCalculationUtilities::EvaluateInPoint(
        this->GetGeometry(), rVariable, rShapeFunction, Step);
}

template <unsigned int TDim, unsigned int TNumNodes>
array_1d<double, 3> SegregatedVMSPressureElement<TDim, TNumNodes>::EvaluateInPoint(
    const Variable<array_1d<double, 3>>& rVariable, const Vector& rShapeFunction, const int Step) const
{
    return RansCalculationUtilities::EvaluateInPoint(
        this->GetGeometry(), rVariable, rShapeFunction, Step);
}

template <unsigned int TDim, unsigned int TNumNodes>
BoundedVector<double, TNumNodes> SegregatedVMSPressureElement<TDim, TNumNodes>::GetConvectionOperator(
    const array_1d<double, 3>& rVector, const Matrix& rShapeDerivatives) const
{
    return SegregatedVMSElementUtilities::GetConvectionOperator<TDim, TNumNodes>(
        rVector, rShapeDerivatives);
}

template <unsigned int TDim, unsigned int TNumNodes>
array_1d<double, 3> SegregatedVMSPressureElement<TDim, TNumNodes>::GetVelocity(
    const Vector& rShapeFunction, const int Step)
{
    const array_1d<double, 3>& r_velocity =
        this->EvaluateInPoint(VELOCITY, rShapeFunction, Step);
    const array_1d<double, 3>& r_mesh_velocity =
        this->EvaluateInPoint(MESH_VELOCITY, rShapeFunction, Step);
    return (r_velocity - r_mesh_velocity);
}

template <unsigned int TDim, unsigned int TNumNodes>
void SegregatedVMSPressureElement<TDim, TNumNodes>::GetFirstDerivativesVector(Vector& Values, int Step)
{
    this->GetValuesVector(Values, Step);
}

template <unsigned int TDim, unsigned int TNumNodes>
void SegregatedVMSPressureElement<TDim, TNumNodes>::EquationIdVector(
    EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
{
    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize);

    for (IndexType i_node = 0; i_node < TNumNodes; ++i_node)
    {
        rResult[i_node] = this->GetGeometry()[i_node].GetDof(PRESSURE).EquationId();
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void SegregatedVMSPressureElement<TDim, TNumNodes>::GetDofList(DofsVectorType& rElementalDofList,
                                                               ProcessInfo& rCurrentProcessInfo)
{
    if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    for (IndexType i_node = 0; i_node < TNumNodes; ++i_node)
    {
        rElementalDofList[i_node] = this->GetGeometry()[i_node].pGetDof(PRESSURE);
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void SegregatedVMSPressureElement<TDim, TNumNodes>::GetValuesVector(Vector& Values, int Step)
{
    if (Values.size() != LocalSize)
        Values.resize(LocalSize, false);

    for (IndexType i_node = 0; i_node < TNumNodes; ++i_node)
    {
        Values[i_node] =
            this->GetGeometry()[i_node].FastGetSolutionStepValue(PRESSURE, Step);
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void SegregatedVMSPressureElement<TDim, TNumNodes>::GetSecondDerivativesVector(Vector& Values, int Step)
{
    if (Values.size() != LocalSize)
        Values.resize(LocalSize, false);

    Values = ZeroVector(LocalSize);
}

template <unsigned int TDim, unsigned int TNumNodes>
int SegregatedVMSPressureElement<TDim, TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo)
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

        KRATOS_CHECK_DOF_IN_NODE(PRESSURE, r_node);
    }

    return value;

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void SegregatedVMSPressureElement<TDim, TNumNodes>::CalculateRightHandSide(
    VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    constexpr IndexType velocity_size = TDim * TNumNodes;

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

    const double dynamic_tau = rCurrentProcessInfo[DYNAMIC_TAU];
    const double delta_time = rCurrentProcessInfo[DELTA_TIME];
    const double element_length = this->GetGeometry().Length();

    const BoundedVector<double, velocity_size>& r_nodal_velocity =
        SegregatedVMSElementUtilities::GetValuesVector<array_1d<double, 3>, velocity_size>(
            this->GetGeometry(), VELOCITY);

    const BoundedVector<double, LocalSize>& r_nodal_pressure =
        SegregatedVMSElementUtilities::GetValuesVector<double, LocalSize>(
            this->GetGeometry(), PRESSURE);

    for (IndexType g = 0; g < num_gauss_points; ++g)
    {
        const Matrix& r_shape_derivatives = shape_derivatives[g];
        const Vector gauss_shape_functions = row(shape_functions, g);
        const double weight = gauss_weights[g];

        const double density = this->EvaluateInPoint(DENSITY, gauss_shape_functions);
        const double dynamic_viscosity = this->EvaluateInPoint(VISCOSITY, gauss_shape_functions) * density;
        const array_1d<double, 3>& r_velocity = this->GetVelocity(gauss_shape_functions);
        const array_1d<double, 3>& r_body_force = this->EvaluateInPoint(BODY_FORCE, gauss_shape_functions) * density;
        const double velocity_magnitude = norm_2(r_velocity);

        const BoundedVector<double, TNumNodes>& r_convective_velocity = this->GetConvectionOperator(r_velocity, r_shape_derivatives);
        const BoundedVector<double, TNumNodes>& r_convective_body_force = this->GetConvectionOperator(r_body_force, r_shape_derivatives);

        const double tau_one = SegregatedVMSElementUtilities::CalculateTauOne(velocity_magnitude, element_length, density, dynamic_viscosity, dynamic_tau, delta_time);

        BoundedMatrix<double, TDim, TDim> velocity_gradient;
        RansCalculationUtilities::CalculateGradient<TDim>(velocity_gradient, this->GetGeometry(), VELOCITY, r_shape_derivatives);
        const double divergence = RansCalculationUtilities::CalculateMatrixTrace<TDim>(velocity_gradient);

        array_1d<double, 3> pressure_gradient;
        RansCalculationUtilities::CalculateGradient(pressure_gradient, this->GetGeometry(), PRESSURE, r_shape_derivatives);
        const BoundedVector<double, TNumNodes>& r_convective_pressure = this->GetConvectionOperator(pressure_gradient, r_shape_derivatives);

        array_1d<double, 3> convective_velocity = ZeroVector(3);
        for (IndexType i = 0; i < TDim; ++i)
        {
            for (IndexType j = 0; j < TDim; ++j)
            {
                convective_velocity[i] += r_velocity[j] * velocity_gradient(i, j);
            }
        }
        const BoundedVector<double, TNumNodes>& r_convective_convective_velocity = this->GetConvectionOperator(convective_velocity, r_shape_derivatives);

        const double coeff_1 = tau_one * weight;
        const double coeff_2 = divergence * weight;
        const double coeff_3 = tau_one * weight * density;

        double value = 0.0;
        for (IndexType a = 0; a < TNumNodes; ++a)
        {
            // adding body force contribution q_f
            value =  coeff_1 * r_convective_body_force[a];

            // adding -D*u
            value -= coeff_2 * gauss_shape_functions[a];

            // adding -Q_u*u
            value -= coeff_3 * r_convective_convective_velocity[a];

            // adding -Q_p*p
            value -= coeff_1 * r_convective_pressure[a];

            rRightHandSideVector[a] += value;
        }
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void SegregatedVMSPressureElement<TDim, TNumNodes>::CalculateLocalVelocityContribution(
    MatrixType& rDampingMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    constexpr IndexType velocity_size = TDim * TNumNodes;

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

    BoundedMatrix<double, TNumNodes, velocity_size> q_u;

    for (IndexType g = 0; g < num_gauss_points; ++g)
    {
        const Matrix& r_shape_derivatives = shape_derivatives[g];
        const Vector& gauss_shape_functions = row(shape_functions, g);
        const double weight = gauss_weights[g];

        const double density = this->EvaluateInPoint(DENSITY, gauss_shape_functions);
        const double dynamic_viscosity = this->EvaluateInPoint(VISCOSITY, gauss_shape_functions) * density;
        const array_1d<double, 3>& r_velocity = this->GetVelocity(gauss_shape_functions);
        const double velocity_magnitude = norm_2(r_velocity);
        const BoundedVector<double, TNumNodes>& r_convective_velocity = this->GetConvectionOperator(r_velocity, r_shape_derivatives);

        const double tau_one = SegregatedVMSElementUtilities::CalculateTauOne(velocity_magnitude, element_length, density, dynamic_viscosity, dynamic_tau, delta_time);
        const double tau_two = SegregatedVMSElementUtilities::CalculateTauTwo(velocity_magnitude, element_length, density, dynamic_viscosity);

        const BoundedMatrix<double, velocity_size, velocity_size>& r_uu =
            SegregatedVMSElementUtilities::CalculateVelocityMatrixGaussPointContributions<TDim, TNumNodes>(
                density, dynamic_viscosity, tau_one, tau_two, r_convective_velocity,
                gauss_shape_functions, r_shape_derivatives, weight);

        BoundedMatrix<double, velocity_size, LocalSize> a_inv_g_sp;
        BoundedMatrix<double, LocalSize, velocity_size> pre_multiplier;

        const double coeff_1 = weight * -1.0;
        const double coeff_2 = density * weight * tau_one ;
        const double coeff_3 = tau_one * weight;

        double value = 0.0;
        for (IndexType a = 0; a < TNumNodes; ++a)
        {
            IndexType row_index = a * TDim;
            for (IndexType b = 0; b < TNumNodes; ++b)
            {
                double dna_dnb = 0.0;
                IndexType col_index = b * TDim;
                for (IndexType i = 0; i < TDim; ++i)
                {
                    // calculating g_sp matrix elements
                    value = coeff_1 * r_shape_derivatives(a, i) * gauss_shape_functions[b];
                    value += coeff_2 * r_convective_velocity[a] *
                             r_shape_derivatives(b, i);
                    // calculating r_uu^-1 * q_sp matrix
                    a_inv_g_sp(row_index + i, b) = value / r_uu(row_index + i, row_index + i);

                    // calculating u' - q_u matrix elements
                    value = r_shape_derivatives(a, i) * gauss_shape_functions[b] * weight;
                    value -= coeff_2 * r_convective_velocity[b] * r_shape_derivatives(a, i);
                    pre_multiplier(a, col_index + i) = value;

                    // calculating q_p matrix elements
                    dna_dnb += r_shape_derivatives(a, i) * r_shape_derivatives(b, i);
                }

                // adding q_p matrix element
                rDampingMatrix(a, b) += coeff_3 * dna_dnb;
            }
        }

        noalias(rDampingMatrix) += prod(pre_multiplier, a_inv_g_sp);
    }

    KRATOS_CATCH("");
}

// template instantiations
template class SegregatedVMSPressureElement<2>;
template class SegregatedVMSPressureElement<3>;

} // namespace Kratos.
