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
#include "includes/checks.h"
#include "includes/variables.h"
#include "includes/cfd_variables.h"

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
        const double dynamic_viscosity =
            this->EvaluateInPoint(VISCOSITY, gauss_shape_functions) * density;
        const array_1d<double, 3>& r_velocity = this->GetVelocity(gauss_shape_functions);
        const double velocity_magnitude = norm_2(r_velocity);
        const double tau_one = SegregatedVMSElementUtilities::CalculateTauOne(
            velocity_magnitude, element_length, density, dynamic_viscosity,
            dynamic_tau, delta_time);
        const array_1d<double, 3>& r_body_force =
            this->EvaluateInPoint(BODY_FORCE, gauss_shape_functions) * density;
        const BoundedVector<double, TNumNodes>& r_convective_body_force =
            this->GetConvectionOperator(r_body_force, r_shape_derivatives);
        const BoundedVector<double, TNumNodes>& r_convective_velocity =
            this->GetConvectionOperator(r_velocity, r_shape_derivatives);

        const double coeff_1 = tau_one * weight;
        const double coeff_2 = tau_one * density * weight;

        // calculating pressure velocity matrix
        BoundedMatrix<double, LocalSize, velocity_size> pu_matrix;

        // calculating Q_p matrix
        BoundedMatrix<double, LocalSize, LocalSize> q_p;

        double value;
        for (IndexType a = 0; a < TNumNodes; ++a)
        {
            for (IndexType b = 0; b < TNumNodes; ++b)
            {
                IndexType col_index = b * TDim;
                double dna_dnb = 0.0;
                for (IndexType d = 0; d < TDim; ++d)
                {
                    dna_dnb += r_shape_derivatives(a, d) * r_shape_derivatives(b, d);

                    value = weight * gauss_shape_functions[a] * r_shape_derivatives(b, d);
                    value += coeff_2 * r_shape_derivatives(a, d) *
                             r_convective_velocity[b];
                    pu_matrix(a, col_index + d) = value;
                }

                q_p(a, b) = coeff_1 * dna_dnb;
            }

            // adding Q_f
            rRightHandSideVector[a] += coeff_1 * r_convective_body_force[a];
        }

        noalias(rRightHandSideVector) -= prod(pu_matrix, r_nodal_velocity);
        noalias(rRightHandSideVector) -= prod(q_p, r_nodal_pressure);
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
        const double dynamic_viscosity =
            this->EvaluateInPoint(VISCOSITY, gauss_shape_functions) * density;
        const array_1d<double, 3>& r_velocity = this->GetVelocity(gauss_shape_functions);
        const double velocity_magnitude = norm_2(r_velocity);
        const double tau_one = SegregatedVMSElementUtilities::CalculateTauOne(
            velocity_magnitude, element_length, density, dynamic_viscosity,
            dynamic_tau, delta_time);
        const double tau_two = SegregatedVMSElementUtilities::CalculateTauTwo(
            velocity_magnitude, element_length, density, dynamic_viscosity);
        const BoundedVector<double, TNumNodes>& r_convective_velocity =
            this->GetConvectionOperator(r_velocity, r_shape_derivatives);

        const double coeff_1 = weight * tau_one;
        const double coeff_2 = weight * tau_one * density;

        for (IndexType a = 0; a < TNumNodes; ++a)
        {
            for (IndexType b = 0; b < TNumNodes; ++b)
            {
                double dna_dnb = 0.0;
                IndexType col_index = b * TDim;

                for (IndexType d = 0; d < TDim; ++d)
                {
                    dna_dnb += r_shape_derivatives(a, d) * r_shape_derivatives(b, d);
                    q_u(a, col_index + d) = (gauss_shape_functions[b] * weight -
                                             coeff_2 * r_convective_velocity[b]) *
                                            r_shape_derivatives(a, d);
                }

                // calculating Q_p
                rDampingMatrix(a, b) += coeff_1 * dna_dnb;
            }
        }

        // Calculating G + S_p matrix
        const BoundedMatrix<double, velocity_size, LocalSize>& r_up_matrix =
            SegregatedVMSElementUtilities::CalculateVelocityPressureMatrixGaussPointContributions<TDim, TNumNodes>(
                density, tau_one, r_convective_velocity, gauss_shape_functions,
                r_shape_derivatives, weight);

        // Calculating C + K + S_u + H_u matrix
        const BoundedMatrix<double, velocity_size, velocity_size>& r_velocity_matrix =
            SegregatedVMSElementUtilities::CalculateVelocityMatrixGaussPointContributions<TDim, TNumNodes>(
                density, dynamic_viscosity, tau_one, tau_two, r_convective_velocity,
                gauss_shape_functions, r_shape_derivatives, weight);

        // Calculating approximated (C + K + S_u + H_u)^(-1) * (G + S_p)
        BoundedMatrix<double, velocity_size, LocalSize> approximated_matrix;
        for (IndexType a = 0; a < velocity_size; ++a)
        {
            // calculating row sum
            double value = 0.0;
            for (IndexType b = 0; b < velocity_size; ++b)
                value += r_velocity_matrix(a, b);

            // assuming row sum is in the diagonal
            for (IndexType b = 0; b < LocalSize; ++b)
                approximated_matrix(a, b) = r_up_matrix(a, b) / value;
        }

        noalias(rDampingMatrix) += prod(q_u, approximated_matrix);
    }

    KRATOS_CATCH("");
}

// template instantiations
template class SegregatedVMSPressureElement<2>;
template class SegregatedVMSPressureElement<3>;

} // namespace Kratos.
