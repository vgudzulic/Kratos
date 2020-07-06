// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:     BSD License
//           license: structural_mechanics_application/license.txt
//
//  Main authors: Ignasi de Pouplana
//
//
//

// System includes

// External includes

// Project includes
#include "custom_elements/truss_element_linear_3D2N.hpp"
#include "includes/define.h"
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"

namespace Kratos {
TrussFICElementLinear3D2N::TrussFICElementLinear3D2N(IndexType NewId,
        GeometryType::Pointer pGeometry)
    : TrussElementLinear3D2N(NewId, pGeometry) {}

TrussFICElementLinear3D2N::TrussFICElementLinear3D2N(
    IndexType NewId, GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties)
    : TrussElementLinear3D2N(NewId, pGeometry, pProperties) {}

Element::Pointer
TrussFICElementLinear3D2N::Create(IndexType NewId,
                               NodesArrayType const& rThisNodes,
                               PropertiesType::Pointer pProperties) const
{
    const GeometryType& rGeom = GetGeometry();
    return Kratos::make_intrusive<TrussFICElementLinear3D2N>(
               NewId, rGeom.Create(rThisNodes), pProperties);
}

Element::Pointer
TrussFICElementLinear3D2N::Create(IndexType NewId,
                               GeometryType::Pointer pGeom,
                               PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<TrussFICElementLinear3D2N>(
               NewId, pGeom, pProperties);
}

TrussFICElementLinear3D2N::~TrussFICElementLinear3D2N() {}


void TrussFICElementLinear3D2N::CalculateRightHandSide(
    VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{

    KRATOS_TRY
    rRightHandSideVector = ZeroVector(msLocalSize);

    BoundedVector<double, msLocalSize> internal_forces = ZeroVector(msLocalSize);
    UpdateInternalForces(internal_forces);

    noalias(rRightHandSideVector) -= internal_forces;

    AddPrestressLinear(rRightHandSideVector);

    // add bodyforces
    noalias(rRightHandSideVector) += CalculateBodyForces();
    KRATOS_CATCH("")
}

void TrussFICElementLinear3D2N::AddExplicitContribution(
    const VectorType& rRHSVector,
    const Variable<VectorType>& rRHSVariable,
    const Variable<double >& rDestinationVariable,
    const ProcessInfo& rCurrentProcessInfo
)
{
    KRATOS_TRY;

    auto& r_geom = GetGeometry();

    if (rDestinationVariable == NODAL_MASS) {
        VectorType element_mass_vector(msLocalSize);
        CalculateLumpedMassVector(element_mass_vector);

        for (SizeType i = 0; i < msNumberOfNodes; ++i) {
            double& r_nodal_mass = r_geom[i].GetValue(NODAL_MASS);
            int index = i * msDimension;

            #pragma omp atomic
            r_nodal_mass += element_mass_vector(index);
        }
    }

    KRATOS_CATCH("")
}

void TrussFICElementLinear3D2N::AddExplicitContribution(
    const VectorType& rRHSVector, const Variable<VectorType>& rRHSVariable,
    const Variable<array_1d<double, 3>>& rDestinationVariable,
    const ProcessInfo& rCurrentProcessInfo
)
{
    KRATOS_TRY;

    if (rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == FORCE_RESIDUAL) {

        BoundedVector<double, msLocalSize> damping_residual_contribution = ZeroVector(msLocalSize);
        Vector current_nodal_velocities = ZeroVector(msLocalSize);
        GetFirstDerivativesVector(current_nodal_velocities);
        Matrix damping_matrix;
        ProcessInfo temp_process_information = rCurrentProcessInfo; // cant pass const ProcessInfo
        CalculateDampingMatrix(damping_matrix, temp_process_information);
        // current residual contribution due to damping
        noalias(damping_residual_contribution) = prod(damping_matrix, current_nodal_velocities);

        for (size_t i = 0; i < msNumberOfNodes; ++i) {
            size_t index = msDimension * i;
            array_1d<double, 3>& r_force_residual = GetGeometry()[i].FastGetSolutionStepValue(FORCE_RESIDUAL);
            for (size_t j = 0; j < msDimension; ++j) {
                #pragma omp atomic
                r_force_residual[j] += rRHSVector[index + j] - damping_residual_contribution[index + j];
            }
        }
    }

    KRATOS_CATCH("")
}


void TrussFICElementLinear3D2N::CalculateDampingMatrixWithLumpedMass(
    MatrixType& rDampingMatrix,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    unsigned int number_of_nodes = GetGeometry().size();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    // Resizing as needed the LHS
    unsigned int mat_size = number_of_nodes * dimension;

    if ( rDampingMatrix.size1() != mat_size )
        rDampingMatrix.resize( mat_size, mat_size, false );

    noalias( rDampingMatrix ) = ZeroMatrix( mat_size, mat_size );

    // 1.-Get Damping Coeffitients (RAYLEIGH_ALPHA, RAYLEIGH_BETA)
    double alpha = 0.0;
    if( GetProperties().Has(RAYLEIGH_ALPHA) )
        alpha = GetProperties()[RAYLEIGH_ALPHA];
    else if( rCurrentProcessInfo.Has(RAYLEIGH_ALPHA) )
        alpha = rCurrentProcessInfo[RAYLEIGH_ALPHA];

    double beta  = 0.0;
    if( GetProperties().Has(RAYLEIGH_BETA) )
        beta = GetProperties()[RAYLEIGH_BETA];
    else if( rCurrentProcessInfo.Has(RAYLEIGH_BETA) )
        beta = rCurrentProcessInfo[RAYLEIGH_BETA];

    // Compose the Damping Matrix:
    // Rayleigh Damping Matrix: alpha*M + beta*K

    // 2.-Calculate mass matrix:
    if (alpha > std::numeric_limits<double>::epsilon()) {
        VectorType temp_vector(mat_size);
        CalculateLumpedMassVector(temp_vector);
        for (IndexType i = 0; i < mat_size; ++i)
            rDampingMatrix(i, i) += alpha * temp_vector[i];
    }

    // 3.-Calculate StiffnessMatrix:
    if (beta > std::numeric_limits<double>::epsilon()) {
        MatrixType stiffness_matrix( mat_size, mat_size );
        VectorType residual_vector( mat_size );

        this->CalculateAll(stiffness_matrix, residual_vector, rCurrentProcessInfo, true, false);

        noalias( rDampingMatrix ) += beta  * stiffness_matrix;
    }

    KRATOS_CATCH( "" )
}

void TrussFICElementLinear3D2N::UpdateInternalForces(BoundedVector<double,msLocalSize>& rInternalForces)
{
    KRATOS_TRY;

    Vector temp_internal_stresses = ZeroVector(msLocalSize);
    ProcessInfo temp_process_information;
    ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),temp_process_information);

    Vector temp_strain = ZeroVector(1);
    Vector temp_stress = ZeroVector(1);
    temp_strain[0] = CalculateLinearStrain();
    Values.SetStrainVector(temp_strain);
    Values.SetStressVector(temp_stress);
    mpConstitutiveLaw->CalculateMaterialResponse(Values,ConstitutiveLaw::StressMeasure_PK2);

    temp_internal_stresses[0] = -1.0*temp_stress[0];
    temp_internal_stresses[3] = temp_stress[0];

    rInternalForces = temp_internal_stresses*GetProperties()[CROSS_AREA];


    BoundedMatrix<double, msLocalSize, msLocalSize> transformation_matrix =
        ZeroMatrix(msLocalSize, msLocalSize);
    CreateTransformationMatrix(transformation_matrix);

    rInternalForces = prod(transformation_matrix, rInternalForces);

    KRATOS_CATCH("");
}


void TrussFICElementLinear3D2N::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, TrussElementLinear3D2N);
    rSerializer.save("mConstitutiveLaw", mpConstitutiveLaw);
}
void TrussFICElementLinear3D2N::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, TrussElementLinear3D2N);
    rSerializer.load("mConstitutiveLaw", mpConstitutiveLaw);
}

} // namespace Kratos.
