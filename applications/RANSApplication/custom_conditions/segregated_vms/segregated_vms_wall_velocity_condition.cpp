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
#include "segregated_vms_wall_velocity_condition.h"

namespace Kratos
{
template <>
void SegregatedVMSWallVelocityCondition<2, 2>::EquationIdVector(EquationIdVectorType& rResult,
                                                                ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int NumNodes = 2;
    const unsigned int LocalSize = 4;
    unsigned int LocalIndex = 0;

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        rResult[LocalIndex++] =
            this->GetGeometry()[iNode].GetDof(VELOCITY_X).EquationId();
        rResult[LocalIndex++] =
            this->GetGeometry()[iNode].GetDof(VELOCITY_Y).EquationId();
    }
}

template <>
void SegregatedVMSWallVelocityCondition<3, 3>::EquationIdVector(EquationIdVectorType& rResult,
                                                                ProcessInfo& rCurrentProcessInfo)
{
    const SizeType NumNodes = 3;
    const SizeType LocalSize = 9;
    unsigned int LocalIndex = 0;

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        rResult[LocalIndex++] =
            this->GetGeometry()[iNode].GetDof(VELOCITY_X).EquationId();
        rResult[LocalIndex++] =
            this->GetGeometry()[iNode].GetDof(VELOCITY_Y).EquationId();
        rResult[LocalIndex++] =
            this->GetGeometry()[iNode].GetDof(VELOCITY_Z).EquationId();
    }
}

template <>
void SegregatedVMSWallVelocityCondition<2, 2>::GetDofList(DofsVectorType& rElementalDofList,
                                                          ProcessInfo& rCurrentProcessInfo)
{
    const SizeType NumNodes = 2;
    const SizeType LocalSize = 4;

    if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    unsigned int LocalIndex = 0;

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_X);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_Y);
    }
}

template <>
void SegregatedVMSWallVelocityCondition<3, 3>::GetDofList(DofsVectorType& rElementalDofList,
                                                          ProcessInfo& rCurrentProcessInfo)
{
    const SizeType NumNodes = 3;
    const SizeType LocalSize = 9;

    if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    unsigned int LocalIndex = 0;

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_X);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_Y);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_Z);
    }
}

template <>
void SegregatedVMSWallVelocityCondition<2, 2>::GetFirstDerivativesVector(Vector& Values, int Step)
{
    const SizeType LocalSize = 4;
    unsigned int LocalIndex = 0;

    if (Values.size() != LocalSize)
        Values.resize(LocalSize, false);

    for (unsigned int iNode = 0; iNode < 2; ++iNode)
    {
        array_1d<double, 3>& r_velocity =
            this->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY, Step);
        for (unsigned int d = 0; d < 2; ++d)
            Values[LocalIndex++] = r_velocity[d];
    }
}

template <>
void SegregatedVMSWallVelocityCondition<3, 3>::GetFirstDerivativesVector(Vector& Values, int Step)
{
    const SizeType LocalSize = 9;
    unsigned int LocalIndex = 0;

    if (Values.size() != LocalSize)
        Values.resize(LocalSize, false);

    for (unsigned int iNode = 0; iNode < 3; ++iNode)
    {
        array_1d<double, 3>& r_velocity =
            this->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY, Step);
        for (unsigned int d = 0; d < 3; ++d)
            Values[LocalIndex++] = r_velocity[d];
    }
}

template <>
void SegregatedVMSWallVelocityCondition<2, 2>::GetSecondDerivativesVector(Vector& Values, int Step)
{
    const SizeType LocalSize = 4;
    unsigned int LocalIndex = 0;

    if (Values.size() != LocalSize)
        Values.resize(LocalSize, false);

    for (unsigned int iNode = 0; iNode < 2; ++iNode)
    {
        array_1d<double, 3>& rVelocity =
            this->GetGeometry()[iNode].FastGetSolutionStepValue(ACCELERATION, Step);
        for (unsigned int d = 0; d < 2; ++d)
            Values[LocalIndex++] = rVelocity[d];
    }
}

template <>
void SegregatedVMSWallVelocityCondition<3, 3>::GetSecondDerivativesVector(Vector& Values, int Step)
{
    const SizeType LocalSize = 9;
    unsigned int LocalIndex = 0;

    if (Values.size() != LocalSize)
        Values.resize(LocalSize, false);

    for (unsigned int iNode = 0; iNode < 3; ++iNode)
    {
        array_1d<double, 3>& rVelocity =
            this->GetGeometry()[iNode].FastGetSolutionStepValue(ACCELERATION, Step);
        for (unsigned int d = 0; d < 3; ++d)
            Values[LocalIndex++] = rVelocity[d];
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
int SegregatedVMSWallVelocityCondition<TDim, TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    int check = BaseType::Check(rCurrentProcessInfo);

    const GeometryType& r_geometry = this->GetGeometry();

    for (IndexType i_node = 0; i_node < TNumNodes; ++i_node)
    {
        const NodeType& r_node = r_geometry[i_node];

        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(MESH_VELOCITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ACCELERATION, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DENSITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(KINEMATIC_VISCOSITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_KINETIC_ENERGY, r_node);

        KRATOS_CHECK_DOF_IN_NODE(VELOCITY_X, r_node);
        KRATOS_CHECK_DOF_IN_NODE(VELOCITY_Y, r_node);
        if (TDim == 3)
            KRATOS_CHECK_DOF_IN_NODE(VELOCITY_Z, r_node);
    }

    return check;

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void SegregatedVMSWallVelocityCondition<TDim, TNumNodes>::Initialize()
{
    KRATOS_TRY;

    if (RansCalculationUtilities::IsWall(*this))
    {
        const ConditionType& r_parent_condition =
            *(this->GetValue(PARENT_CONDITION_POINTER));

        const array_1d<double, 3>& rNormal = r_parent_condition.GetValue(NORMAL);
        KRATOS_ERROR_IF(norm_2(rNormal) == 0.0)
            << "NORMAL must be calculated before using this " << this->Info() << "\n";

        KRATOS_ERROR_IF(r_parent_condition.GetValue(NEIGHBOUR_ELEMENTS).size() == 0)
            << this->Info() << " cannot find parent element\n";

        const double nu = RansCalculationUtilities::EvaluateInParentCenter(
            KINEMATIC_VISCOSITY, r_parent_condition);
        KRATOS_ERROR_IF(nu == 0.0)
            << "KINEMATIC_VISCOSITY is not defined in the parent element of "
            << this->Info() << "\n.";

        mWallHeight =
            RansCalculationUtilities::CalculateWallHeight(r_parent_condition, rNormal);
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
std::string SegregatedVMSWallVelocityCondition<TDim, TNumNodes>::Info() const
{
    std::stringstream buffer;
    buffer << "SegregatedVMSWallVelocityCondition" << TDim << "D";
    return buffer.str();
}

template <unsigned int TDim, unsigned int TNumNodes>
void SegregatedVMSWallVelocityCondition<TDim, TNumNodes>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "SegregatedVMSWallVelocityCondition";
}

template <unsigned int TDim, unsigned int TNumNodes>
void SegregatedVMSWallVelocityCondition<TDim, TNumNodes>::PrintData(std::ostream& rOStream) const
{
}

template <unsigned int TDim, unsigned int TNumNodes>
void SegregatedVMSWallVelocityCondition<TDim, TNumNodes>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    const SizeType BlockSize = TDim;
    const SizeType LocalSize = BlockSize * TNumNodes;

    if (rLeftHandSideMatrix.size1() != LocalSize)
        rLeftHandSideMatrix.resize(LocalSize, LocalSize);

    if (rRightHandSideVector.size() != LocalSize)
        rRightHandSideVector.resize(LocalSize);

    noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize);
    noalias(rRightHandSideVector) = ZeroVector(LocalSize);
}

template <unsigned int TDim, unsigned int TNumNodes>
void SegregatedVMSWallVelocityCondition<TDim, TNumNodes>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{
    const SizeType BlockSize = TDim;
    const SizeType LocalSize = BlockSize * TNumNodes;

    if (rLeftHandSideMatrix.size1() != LocalSize)
        rLeftHandSideMatrix.resize(LocalSize, LocalSize);

    noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize);
}

template <unsigned int TDim, unsigned int TNumNodes>
void SegregatedVMSWallVelocityCondition<TDim, TNumNodes>::CalculateRightHandSide(
    VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    const SizeType BlockSize = TDim;
    const SizeType LocalSize = BlockSize * TNumNodes;

    if (rRightHandSideVector.size() != LocalSize)
        rRightHandSideVector.resize(LocalSize);

    noalias(rRightHandSideVector) = ZeroVector(LocalSize);
}

template <unsigned int TDim, unsigned int TNumNodes>
void SegregatedVMSWallVelocityCondition<TDim, TNumNodes>::CalculateDampingMatrix(
    MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo)
{
    VectorType RHS;
    this->CalculateLocalVelocityContribution(rDampingMatrix, RHS, rCurrentProcessInfo);
}

template <unsigned int TDim, unsigned int TNumNodes>
void SegregatedVMSWallVelocityCondition<TDim, TNumNodes>::CalculateLocalVelocityContribution(
    MatrixType& rDampMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    // Initialize local contributions
    const SizeType LocalSize = TDim * TNumNodes;

    if (rDampMatrix.size1() != LocalSize)
        rDampMatrix.resize(LocalSize, LocalSize);
    if (rRightHandSideVector.size() != LocalSize)
        rRightHandSideVector.resize(LocalSize);

    noalias(rDampMatrix) = ZeroMatrix(LocalSize, LocalSize);
    noalias(rRightHandSideVector) = ZeroVector(LocalSize);

    this->ApplyWallLaw(rDampMatrix, rRightHandSideVector, rCurrentProcessInfo);
}

template <unsigned int TDim, unsigned int TNumNodes>
void SegregatedVMSWallVelocityCondition<TDim, TNumNodes>::ApplyWallLaw(
    MatrixType& rLocalMatrix, VectorType& rLocalVector, ProcessInfo& rCurrentProcessInfo)
{
    if (RansCalculationUtilities::IsWall(*this))
    {
        const double eps = std::numeric_limits<double>::epsilon();

        const ConditionType& r_parent_condition =
            *(this->GetValue(PARENT_CONDITION_POINTER));

        const array_1d<double, 3> wall_cell_center_velocity =
            RansCalculationUtilities::CalculateWallVelocity(r_parent_condition);
        const double wall_cell_center_velocity_magnitude = norm_2(wall_cell_center_velocity);

        const double y_plus_limit = rCurrentProcessInfo[RANS_Y_PLUS_LIMIT];
        double y_plus = 0.0;

        if (wall_cell_center_velocity_magnitude > eps)
        {
            constexpr unsigned int block_size = TDim;

            // calculate cell centered y_plus value
            const double kappa = rCurrentProcessInfo[WALL_VON_KARMAN];
            const double beta = rCurrentProcessInfo[WALL_SMOOTHNESS_BETA];
            const double nu = RansCalculationUtilities::EvaluateInParentCenter(
                KINEMATIC_VISCOSITY, r_parent_condition);

            double u_tau;
            RansCalculationUtilities::CalculateYPlusAndUtau(
                y_plus, u_tau, wall_cell_center_velocity_magnitude, mWallHeight,
                nu, kappa, beta);

            GeometryType& r_geometry = this->GetGeometry();

            MatrixType shape_functions;
            VectorType gauss_weights;
            RansCalculationUtilities::CalculateConditionGeometryData(
                r_geometry, GeometryData::GI_GAUSS_2, gauss_weights, shape_functions);
            const int number_of_gauss_points = gauss_weights.size();

            // In the linear region, force the velocity to be in the log region lowest
            // since k - epsilon is only valid in the log region.
            // In order to avoid issues with stagnation points, tke is also used
            const double c_mu_25 =
                std::pow(rCurrentProcessInfo[TURBULENCE_RANS_C_MU], 0.25);
            const std::function<double(double, double)> linear_region_functional =
                [c_mu_25, y_plus_limit](const double TurbulentKineticEnergy,
                                        const double Velocity) -> double {
                return std::max(c_mu_25 * std::sqrt(std::max(TurbulentKineticEnergy, 0.0)),
                                Velocity / y_plus_limit);
            };

            // log region, apply the u_tau which is calculated based on the cell centered velocity
            const std::function<double(double, double)> log_region_functional =
                [u_tau](const double, const double) -> double { return u_tau; };

            const std::function<double(double, double)> wall_tau_function =
                ((y_plus >= y_plus_limit) ? log_region_functional : linear_region_functional);

            double condition_u_tau = 0.0;

            for (int g = 0; g < number_of_gauss_points; ++g)
            {
                const Vector& gauss_shape_functions = row(shape_functions, g);
                const double weight = gauss_weights[g];

                const array_1d<double, 3>& r_wall_velocity =
                    RansCalculationUtilities::EvaluateInPoint(
                        r_geometry, VELOCITY, gauss_shape_functions);
                const double wall_velocity_magnitude = norm_2(r_wall_velocity);
                const double rho = RansCalculationUtilities::EvaluateInPoint(
                    r_geometry, DENSITY, gauss_shape_functions);

                if (wall_velocity_magnitude > eps)
                {
                    const double tke = RansCalculationUtilities::EvaluateInPoint(
                        r_geometry, TURBULENT_KINETIC_ENERGY, gauss_shape_functions);
                    const double gauss_u_tau =
                        wall_tau_function(tke, wall_velocity_magnitude);

                    condition_u_tau += gauss_u_tau;

                    const double coeff_1 = rho * std::pow(gauss_u_tau, 2) *
                                           weight / wall_velocity_magnitude;

                    for (IndexType a = 0; a < TNumNodes; ++a)
                    {
                        for (IndexType i = 0; i < TDim; ++i)
                        {
                            for (IndexType b = 0; b < TNumNodes; ++b)
                            {
                                rLocalMatrix(a * block_size + i, b * block_size + i) +=
                                    gauss_shape_functions[a] *
                                    gauss_shape_functions[b] * coeff_1;
                            }
                            rLocalVector[a * block_size + i] -=
                                gauss_shape_functions[a] * coeff_1 * r_wall_velocity[i];
                        }
                    }
                }
            }

            condition_u_tau /= static_cast<double>(number_of_gauss_points);
            this->GetValue(PARENT_CONDITION_POINTER)->SetValue(FRICTION_VELOCITY, condition_u_tau);
        }

        this->GetValue(PARENT_CONDITION_POINTER)->SetValue(RANS_Y_PLUS, std::max(y_plus, y_plus_limit));
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void SegregatedVMSWallVelocityCondition<TDim, TNumNodes>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
}

template <unsigned int TDim, unsigned int TNumNodes>
void SegregatedVMSWallVelocityCondition<TDim, TNumNodes>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
}

// template instantiations

template class SegregatedVMSWallVelocityCondition<2, 2>;
template class SegregatedVMSWallVelocityCondition<3, 3>;

} // namespace Kratos.
