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
#include "custom_elements/truss_fic_element_linear_3D2N.hpp"
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

        VectorType element_damping_vector(msLocalSize);
        CalculateLumpedDampingVector(element_damping_vector, rCurrentProcessInfo);

        VectorType element_mass_vector(msLocalSize);
        CalculateLumpedMassVector(element_mass_vector);

        for (SizeType i = 0; i < msNumberOfNodes; ++i) {
            double& r_nodal_damping = r_geom[i].GetValue(NODAL_DISPLACEMENT_DAMPING);
            double& r_nodal_mass = r_geom[i].GetValue(NODAL_MASS);
            int index = i * msDimension;

            #pragma omp atomic
            r_nodal_damping += element_damping_vector[index];

            #pragma omp atomic
            r_nodal_mass += rCurrentProcessInfo[DELTA_TIME]*1.0 +
                            rCurrentProcessInfo[MASS_FACTOR]*element_mass_vector[index]/element_damping_vector[index];


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

        VectorType element_damping_vector(msLocalSize);
        CalculateLumpedDampingVector(element_damping_vector, rCurrentProcessInfo);

        VectorType element_mass_vector(msLocalSize);
        CalculateLumpedMassVector(element_mass_vector);

        Vector current_nodal_displacements = ZeroVector(msLocalSize);
        GetValuesVector(current_nodal_displacements);

        for (size_t i = 0; i < msNumberOfNodes; ++i) {
            size_t index = msDimension * i;
            array_1d<double, 3>& r_force_residual = GetGeometry()[i].FastGetSolutionStepValue(FORCE_RESIDUAL);
            array_1d<double, 3>& r_inertial_residual = GetGeometry()[i].FastGetSolutionStepValue(NODAL_INERTIA);
            for (size_t j = 0; j < msDimension; ++j) {
                #pragma omp atomic
                r_force_residual[j] += rRHSVector[index + j];

                #pragma omp atomic
                r_inertial_residual[j] += rCurrentProcessInfo[MASS_FACTOR]*element_mass_vector[index + j]*current_nodal_displacements[index + j]/element_damping_vector[index + j];
            }
        }
    }

    KRATOS_CATCH("")
}

void TrussFICElementLinear3D2N::CalculateLumpedMassVector(VectorType& rMassVector)
{
    KRATOS_TRY

    // Clear matrix
    if (rMassVector.size() != msLocalSize) {
        rMassVector.resize(msLocalSize, false);
    }

    const double A = GetProperties()[CROSS_AREA];
    const double L = StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this);
    const double rho = GetProperties()[DENSITY];

    const double total_mass = A * L * rho;

    for (int i = 0; i < msNumberOfNodes; ++i) {
        for (int j = 0; j < msDimension; ++j) {
            int index = i * msDimension + j;

            rMassVector[index] = total_mass * 0.50;
        }
    }

    KRATOS_CATCH("")
}

void TrussFICElementLinear3D2N::CalculateLumpedStiffnessVector(VectorType& rStiffnessVector,const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Clear Vector
    if (rStiffnessVector.size() != msLocalSize) {
        rStiffnessVector.resize(msLocalSize, false);
    }

    MatrixType stiffness_matrix( msLocalSize, msLocalSize );
    noalias(stiffness_matrix) = ZeroMatrix(msLocalSize,msLocalSize);
    ProcessInfo temp_process_information = rCurrentProcessInfo;
    noalias(stiffness_matrix) = CreateElementStiffnessMatrix(temp_process_information);
    // TODO: this is a first approximation
    for (IndexType i = 0; i < msLocalSize; ++i)
        rStiffnessVector[i] = stiffness_matrix(i,i);

    KRATOS_CATCH("")
}

void TrussFICElementLinear3D2N::CalculateLumpedDampingVector(
    VectorType& rDampingVector,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    // Clear Vector
    if (rDampingVector.size() != msLocalSize) {
        rDampingVector.resize(msLocalSize, false);
    }
    noalias(rDampingVector) = ZeroVector(msLocalSize);

    // 1.-Get Damping Coefficients (RAYLEIGH_ALPHA, RAYLEIGH_BETA)
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

    // Compose the Damping Vector:
    // Rayleigh Damping Vector: alpha*M + beta*K

    // 2.-Calculate mass Vector:
    if (alpha > std::numeric_limits<double>::epsilon()) {
        VectorType mass_vector(msLocalSize);
        CalculateLumpedMassVector(mass_vector);
        for (IndexType i = 0; i < msLocalSize; ++i)
            rDampingVector[i] += alpha * mass_vector[i];

        KRATOS_WATCH(mass_vector)
    }

    // 3.-Calculate Stiffness Vector:
    if (beta > std::numeric_limits<double>::epsilon()) {
        VectorType stiffness_vector(msLocalSize);
        CalculateLumpedStiffnessVector(stiffness_vector,rCurrentProcessInfo);
        for (IndexType i = 0; i < msLocalSize; ++i)
            rDampingVector[i] += beta * stiffness_vector[i];

        KRATOS_WATCH(stiffness_vector)
    }

    KRATOS_WATCH(rDampingVector)

    KRATOS_CATCH( "" )
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
