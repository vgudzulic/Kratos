// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo
//
//

// System includes
#include <iostream>

// External includes

// Project includes
#include "includes/checks.h"
#include "includes/properties.h"
#include "custom_advanced_constitutive/hyper_elastic_isotropic_kirchhoff_3d.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
//******************************CONSTRUCTOR*******************************************
/***********************************************************************************/
HyperElasticIsotropicOgden3D::HyperElasticIsotropicOgden3D()
    : ConstitutiveLaw()
{
}

//******************************COPY CONSTRUCTOR**************************************
/***********************************************************************************/

HyperElasticIsotropicOgden3D::HyperElasticIsotropicOgden3D(const HyperElasticIsotropicOgden3D& rOther)
    : ConstitutiveLaw(rOther)
{
}

//********************************CLONE***********************************************
/***********************************************************************************/

ConstitutiveLaw::Pointer HyperElasticIsotropicOgden3D::Clone() const
{
    return Kratos::make_shared<HyperElasticIsotropicOgden3D>(*this);
}

//*******************************DESTRUCTOR*******************************************
/***********************************************************************************/

HyperElasticIsotropicOgden3D::~HyperElasticIsotropicOgden3D()
{
};

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicOgden3D::CalculateMaterialResponsePK1 (ConstitutiveLaw::Parameters& rValues)
{
    HyperElasticIsotropicOgden3D::CalculateMaterialResponsePK2(rValues);

    Vector& stress_vector                = rValues.GetStressVector();
    const Matrix& deformation_gradient_f = rValues.GetDeformationGradientF();
    const double determinant_f           = rValues.GetDeterminantF();

    this->TransformStresses(stress_vector, deformation_gradient_f, determinant_f, StressMeasure_PK2, StressMeasure_PK1);
}

/***********************************************************************************/
/***********************************************************************************/

void  HyperElasticIsotropicOgden3D::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicOgden3D::CalculateMaterialResponseKirchhoff (ConstitutiveLaw::Parameters& rValues)
{

}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicOgden3D::CalculateMaterialResponseCauchy (ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY;

    // Get Values to compute the constitutive law:
    Flags& r_flags= rValues.GetOptions();

    const Properties& material_properties = rValues.GetMaterialProperties();















    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicOgden3D::InitializeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicOgden3D::InitializeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicOgden3D::InitializeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicOgden3D::InitializeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicOgden3D::FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicOgden3D::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicOgden3D::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicOgden3D::FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
}

/***********************************************************************************/
/***********************************************************************************/

double& HyperElasticIsotropicOgden3D::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    const Properties& material_properties  = rParameterValues.GetMaterialProperties();
    Vector& strain_vector                  = rParameterValues.GetStrainVector();

    // The material properties
    const double young_modulus = material_properties[YOUNG_MODULUS];
    const double poisson_coefficient = material_properties[POISSON_RATIO];

    // The LAME parameters
    const double lame_lambda = (young_modulus * poisson_coefficient)/((1.0 + poisson_coefficient)*(1.0 - 2.0 * poisson_coefficient));
    const double lame_mu = young_modulus/(2.0 * (1.0 + poisson_coefficient));

    // We compute the right Cauchy-Green tensor (C):
    if (rThisVariable == STRAIN_ENERGY) {
        CalculateGreenLagrangianStrain(rParameterValues,strain_vector);
        const Matrix E_tensor=MathUtils<double>::StrainVectorToTensor(strain_vector);
        const Matrix E_tensor_sq=prod(E_tensor,E_tensor);
        double E_trace = 0.0;
        double E_trace_sq = 0.0;
        for (IndexType i = 0; i < E_tensor.size1();i++) {
            E_trace     += E_tensor (i,i);
            E_trace_sq  += E_tensor_sq(i,i);
        }

        rValue = 0.5 * lame_lambda * E_trace * E_trace  + 0.5 *lame_mu *E_trace_sq ;
    }

    return( rValue );
}

/***********************************************************************************/
/***********************************************************************************/

Vector& HyperElasticIsotropicOgden3D::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    if (rThisVariable == STRAIN ||
        rThisVariable == GREEN_LAGRANGE_STRAIN_VECTOR ||
        rThisVariable == HENCKY_STRAIN_VECTOR ||
        rThisVariable == BIOT_STRAIN_VECTOR ||
        rThisVariable == ALMANSI_STRAIN_VECTOR) {

        // Get Values to compute the constitutive law:
        Flags& r_flags = rParameterValues.GetOptions();

        // Previous flags saved
        const bool flag_strain = r_flags.Is( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN );
        const bool flag_const_tensor = r_flags.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR );
        const bool flag_stress = r_flags.Is( ConstitutiveLaw::COMPUTE_STRESS );

        r_flags.Set( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false );
        r_flags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false );
        r_flags.Set( ConstitutiveLaw::COMPUTE_STRESS, false );

        // We compute the strain
        if (rThisVariable == STRAIN) {
            HyperElasticIsotropicOgden3D::CalculateMaterialResponse(rParameterValues, this->GetStressMeasure());
        } else if (rThisVariable == GREEN_LAGRANGE_STRAIN_VECTOR) {
            HyperElasticIsotropicOgden3D::CalculateMaterialResponsePK2(rParameterValues);
        } else if (rThisVariable == ALMANSI_STRAIN_VECTOR) {
            HyperElasticIsotropicOgden3D::CalculateMaterialResponseKirchhoff(rParameterValues);
        } else if (rThisVariable == HENCKY_STRAIN_VECTOR) {
            const Matrix& deformation_gradient_f = rParameterValues.GetDeformationGradientF();
            const Matrix C_tensor = prod(trans( deformation_gradient_f), deformation_gradient_f);
            Vector& r_strain_vector = rParameterValues.GetStrainVector();
            ConstitutiveLawUtilities<VoigtSize>::CalculateHenckyStrain(C_tensor, r_strain_vector);
        } else if (rThisVariable == BIOT_STRAIN_VECTOR) {
            const Matrix& deformation_gradient_f = rParameterValues.GetDeformationGradientF();
            const Matrix C_tensor = prod(trans( deformation_gradient_f), deformation_gradient_f);
            Vector& r_strain_vector = rParameterValues.GetStrainVector();
            ConstitutiveLawUtilities<VoigtSize>::CalculateBiotStrain(C_tensor, r_strain_vector);
        }

        rValue = rParameterValues.GetStrainVector();

        // Previous flags restored
        r_flags.Set( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain );
        r_flags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor );
        r_flags.Set( ConstitutiveLaw::COMPUTE_STRESS, flag_stress );
    } else if (rThisVariable == STRESSES ||
        rThisVariable == CAUCHY_STRESS_VECTOR ||
        rThisVariable == KIRCHHOFF_STRESS_VECTOR ||
        rThisVariable == PK2_STRESS_VECTOR) {

        // Get Values to compute the constitutive law:
        Flags& r_flags = rParameterValues.GetOptions();

        // Previous flags saved
        const bool flag_strain = r_flags.Is( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN );
        const bool flag_const_tensor = r_flags.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR );
        const bool flag_stress = r_flags.Is( ConstitutiveLaw::COMPUTE_STRESS );

        r_flags.Set( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false );
        r_flags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true );
        r_flags.Set( ConstitutiveLaw::COMPUTE_STRESS, true );

        // We compute the stress
        if (rThisVariable == STRESSES) {
            HyperElasticIsotropicOgden3D::CalculateMaterialResponse(rParameterValues, this->GetStressMeasure());
        } if (rThisVariable == KIRCHHOFF_STRESS_VECTOR) {
            HyperElasticIsotropicOgden3D::CalculateMaterialResponseKirchhoff(rParameterValues);
        } if (rThisVariable == CAUCHY_STRESS_VECTOR) {
            HyperElasticIsotropicOgden3D::CalculateMaterialResponseCauchy(rParameterValues);
        } if (rThisVariable == PK2_STRESS_VECTOR) {
            HyperElasticIsotropicOgden3D::CalculateMaterialResponsePK2(rParameterValues);
        }

        rValue = rParameterValues.GetStressVector();

        // Previous flags restored
        r_flags.Set( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain );
        r_flags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor );
        r_flags.Set( ConstitutiveLaw::COMPUTE_STRESS, flag_stress );
    }

    return( rValue );
}

/***********************************************************************************/
/***********************************************************************************/

Matrix& HyperElasticIsotropicOgden3D::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    if (rThisVariable == CONSTITUTIVE_MATRIX ||
        rThisVariable == CONSTITUTIVE_MATRIX_PK2 ||
        rThisVariable == CONSTITUTIVE_MATRIX_KIRCHHOFF) {
        // Get Values to compute the constitutive law:
        Flags& r_flags = rParameterValues.GetOptions();

        // Previous flags saved
        const bool flag_strain = r_flags.Is( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN );
        const bool flag_const_tensor = r_flags.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR );
        const bool flag_stress = r_flags.Is( ConstitutiveLaw::COMPUTE_STRESS );

        r_flags.Set( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false );
        r_flags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true );
        r_flags.Set( ConstitutiveLaw::COMPUTE_STRESS, false );

        // We compute the constitutive matrix
        if (rThisVariable == CONSTITUTIVE_MATRIX) {
            HyperElasticIsotropicOgden3D::CalculateMaterialResponse(rParameterValues, this->GetStressMeasure());
        } else if (rThisVariable == CONSTITUTIVE_MATRIX_PK2) {
            HyperElasticIsotropicOgden3D::CalculateMaterialResponsePK2(rParameterValues);
        } else if (rThisVariable == CONSTITUTIVE_MATRIX_KIRCHHOFF) {
            HyperElasticIsotropicOgden3D::CalculateMaterialResponsePK2(rParameterValues);
        }

        rValue = rParameterValues.GetConstitutiveMatrix();

        // Previous flags restored
        r_flags.Set( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain );
        r_flags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor );
        r_flags.Set( ConstitutiveLaw::COMPUTE_STRESS, flag_stress );
    }

    return( rValue );
}

//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
/***********************************************************************************/

void HyperElasticIsotropicOgden3D::GetLawFeatures(Features& rFeatures)
{
    //Set the type of law
    rFeatures.mOptions.Set( THREE_DIMENSIONAL_LAW );
    rFeatures.mOptions.Set( FINITE_STRAINS );
    rFeatures.mOptions.Set( ISOTROPIC );

    //Set strain measure required by the consitutive law
    rFeatures.mStrainMeasures.push_back(StrainMeasure_GreenLagrange);
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

    //Set the strain size
    rFeatures.mStrainSize = GetStrainSize();

    //Set the spacedimension
    rFeatures.mSpaceDimension = WorkingSpaceDimension();
}

/***********************************************************************************/
/***********************************************************************************/

int HyperElasticIsotropicOgden3D::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_CHECK_VARIABLE_KEY(YOUNG_MODULUS);
    KRATOS_ERROR_IF(rMaterialProperties[YOUNG_MODULUS] <= 0.0) << "YOUNG_MODULUS is null or negative." << std::endl;

    KRATOS_CHECK_VARIABLE_KEY(POISSON_RATIO);
    const double tolerance = 1.0e-12;
    const double nu_upper_bound = 0.5;
    const double nu_lower_bound = -1.0;
    const double nu = rMaterialProperties[POISSON_RATIO];
    KRATOS_ERROR_IF((nu_upper_bound - nu) < tolerance) << "POISSON_RATIO is above the upper bound 0.5." << std::endl;
    KRATOS_ERROR_IF((nu - nu_lower_bound) < tolerance) << "POISSON_RATIO is below the lower bound -1.0." << std::endl;

    KRATOS_CHECK_VARIABLE_KEY(DENSITY);
    KRATOS_ERROR_IF(rMaterialProperties[DENSITY] < 0.0) << "DENSITY is negative." << std::endl;

    return 0;
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicOgden3D::CalculateConstitutiveMatrixPK2(
    Matrix& rConstitutiveMatrix,
    const double YoungModulus,
    const double PoissonCoefficient
    )
{
    rConstitutiveMatrix.clear();
    rConstitutiveMatrix = ZeroMatrix(6,6);
    const double c1 = YoungModulus / (( 1.00 + PoissonCoefficient ) * ( 1 - 2 * PoissonCoefficient ) );
    const double c2 = c1 * ( 1 - PoissonCoefficient );
    const double c3 = c1 * PoissonCoefficient;
    const double c4 = c1 * 0.5 * ( 1 - 2 * PoissonCoefficient );

    rConstitutiveMatrix( 0, 0 ) = c2;
    rConstitutiveMatrix( 0, 1 ) = c3;
    rConstitutiveMatrix( 0, 2 ) = c3;
    rConstitutiveMatrix( 1, 0 ) = c3;
    rConstitutiveMatrix( 1, 1 ) = c2;
    rConstitutiveMatrix( 1, 2 ) = c3;
    rConstitutiveMatrix( 2, 0 ) = c3;
    rConstitutiveMatrix( 2, 1 ) = c3;
    rConstitutiveMatrix( 2, 2 ) = c2;
    rConstitutiveMatrix( 3, 3 ) = c4;
    rConstitutiveMatrix( 4, 4 ) = c4;
    rConstitutiveMatrix( 5, 5 ) = c4;
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicOgden3D::CalculateKirchhoffStress(
    const Vector& rStrainVector,
    Vector& rStressVector,
    const Matrix& rDeformationGradientF,
    const double YoungModulus,
    const double PoissonCoefficient
    )
{
    CalculatePK2Stress( rStrainVector, rStressVector, YoungModulus, PoissonCoefficient );
    Matrix stress_matrix = MathUtils<double>::StressVectorToTensor( rStressVector );
    ContraVariantPushForward (stress_matrix,rDeformationGradientF); //Kirchhoff
    rStressVector = MathUtils<double>::StressTensorToVector( stress_matrix, rStressVector.size() );
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicOgden3D::CalculatePK2Stress(
    const Vector& rStrainVector,
    Vector& rStressVector,
    const double YoungModulus,
    const double PoissonCoefficient
    )
{
    const double lame_lambda = (YoungModulus * PoissonCoefficient)/((1.0 + PoissonCoefficient)*(1.0 - 2.0 * PoissonCoefficient));
    const double lame_mu = YoungModulus/(2.0 * (1.0 + PoissonCoefficient));
    const Matrix E_tensor=MathUtils<double>::StrainVectorToTensor(rStrainVector);
    double E_trace = 0.0;
    for (unsigned int i = 0; i < E_tensor.size1();i++) {
      E_trace += E_tensor (i,i);
    }
    const SizeType dimension = WorkingSpaceDimension();
    Matrix stress_matrix = lame_lambda*E_trace*IdentityMatrix(dimension) + 2.0 * lame_mu * E_tensor;
    rStressVector = MathUtils<double>::StressTensorToVector( stress_matrix, rStressVector.size() );

//     // Other possibility
//     SizeType size_system = GetStrainSize();
//     Matrix constitutive_matrix_pk2(size_system,size_system);
//     CalculateConstitutiveMatrixPK2(constitutive_matrix_pk2, YoungModulus, PoissonCoefficient);
//     noalias(rStressVector) = prod(constitutive_matrix_pk2,rStrainVector);
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicOgden3D::CalculateConstitutiveMatrixKirchhoff(
    Matrix& ConstitutiveMatrix,
    const Matrix& DeformationGradientF,
    const double YoungModulus,
    const double PoissonCoefficient
    )
{
    ConstitutiveMatrix.clear();

    this->CalculateConstitutiveMatrixPK2(ConstitutiveMatrix, YoungModulus, PoissonCoefficient);
    PushForwardConstitutiveMatrix (ConstitutiveMatrix,DeformationGradientF );
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicOgden3D::CalculateGreenLagrangianStrain(
    ConstitutiveLaw::Parameters& rValues,
    Vector& rStrainVector
    )
{
    // 1.-Compute total deformation gradient
    const Matrix& F = rValues.GetDeformationGradientF();

    // 2.-Compute e = 0.5*(inv(C) - I)
    const Matrix C_tensor = prod(trans(F),F);
    ConstitutiveLawUtilities<VoigtSize>::CalculateGreenLagrangianStrain(C_tensor, rStrainVector);
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicOgden3D::CalculateAlmansiStrain(
    ConstitutiveLaw::Parameters& rValues,
    Vector& rStrainVector
    )
{
    // 1.-Compute total deformation gradient
    const Matrix& F = rValues.GetDeformationGradientF();

    // 2.-COmpute e = 0.5*(1-inv(B))
    const Matrix B_tensor = prod(F,trans(F));
    ConstitutiveLawUtilities<VoigtSize>::CalculateAlmansiStrain(B_tensor, rStrainVector);
}

} // Namespace Kratos
