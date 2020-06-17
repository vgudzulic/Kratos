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
#include "custom_advanced_constitutive/hyper_elastic_isotropic_ogden_3d.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/tangent_operator_calculator_utility.h"

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

void  HyperElasticIsotropicOgden3D::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY;

    HyperElasticIsotropicOgden3D::CalculateMaterialResponsePK2(rValues);

    Vector& stress_vector                = rValues.GetStressVector();
    const Matrix& deformation_gradient_f = rValues.GetDeformationGradientF();
    const double determinant_f           = rValues.GetDeterminantF();

    this->TransformStresses(stress_vector, deformation_gradient_f, determinant_f, StressMeasure_PK2, StressMeasure_Cauchy);

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicOgden3D::CalculateMaterialResponseKirchhoff (ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY;

    HyperElasticIsotropicOgden3D::CalculateMaterialResponsePK2(rValues);

    Vector& stress_vector                = rValues.GetStressVector();
    const Matrix& deformation_gradient_f = rValues.GetDeformationGradientF();
    const double determinant_f           = rValues.GetDeterminantF();

    this->TransformStresses(stress_vector, deformation_gradient_f, determinant_f, StressMeasure_PK2, StressMeasure_Kirchhoff);

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicOgden3D::CalculateMaterialResponsePK2 (ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY;

    // Get Values to compute the constitutive law:
    Flags& r_flags= rValues.GetOptions();

    const Properties& r_material_properties = rValues.GetMaterialProperties();
    Vector& r_strain_vector                   = rValues.GetStrainVector();

    if (r_flags.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        this->CalculateGreenLagrangianStrain(rValues, r_strain_vector);
    }

    if (r_flags.Is( ConstitutiveLaw::COMPUTE_STRESS)) {
        Vector& r_stress_vector = rValues.GetStressVector();
        const Vector& r_ogden_parameters = r_material_properties[OGDEN_HYPERELASTIC_PARAMETERS];

        const double mu_1    = r_ogden_parameters(0);
        const double mu_2    = r_ogden_parameters(1);
        const double mu_3    = r_ogden_parameters(2);
        const double alpha_1 = r_ogden_parameters(3);
        const double alpha_2 = r_ogden_parameters(4);
        const double alpha_3 = r_ogden_parameters(5);

        const Matrix& rFDeformationGradient = rValues.GetDeformationGradientF();

        // We compute Right Cauchy tensor
        const Matrix& r_C = prod(trans(rFDeformationGradient), rFDeformationGradient);

        // Decompose matrix C
        BoundedMatrix<double, Dimension, Dimension> eigen_vector_matrix, eigen_values_matrix, deviatoric_eigen_values;
        MathUtils<double>::GaussSeidelEigenSystem(r_C, eigen_vector_matrix, eigen_values_matrix, 1.0e-16, 200);

        noalias(deviatoric_eigen_values) = eigen_values_matrix * 
            std::pow(eigen_values_matrix(0, 0) * eigen_values_matrix(1, 1) * eigen_values_matrix(2, 2), -1.0 / 3.0);

        std::vector<Vector> eigen_vectors_container(Dimension);
        std::vector<Matrix> M_container(Dimension);

        // We fill the eigen_vectors_container
        Vector auxiliar_vector = ZeroVector(Dimension);
        for (IndexType i = 0; i < Dimension; ++i) {
            for (IndexType j = 0; j < Dimension; ++j) {
                auxiliar_vector[j] = eigen_vector_matrix(j, i);
            }
            noalias(eigen_vectors_container[i]) = auxiliar_vector;
        }

        // We fill the M_container
        Matrix auxiliar_matrix(Dimension, Dimension);
        for (IndexType i = 0; i < Dimension; ++i) {
            noalias(auxiliar_matrix) = std::pow(eigen_values_matrix(i, i), -2) * outer_prod(eigen_vectors_container[i], eigen_vectors_container[i]);
            noalias(M_container[i])  = auxiliar_matrix;
        }

        // We fill the beta parameter
        array_1d<double, Dimension> beta;
        for (IndexType i = 0; i < Dimension; ++i) {
            beta(i) = mu_1 * (std::pow(deviatoric_eigen_values(i, i), alpha_1) - (std::pow(deviatoric_eigen_values(0, 0), alpha_1) + std::pow(deviatoric_eigen_values(1, 1), alpha_1) + std::pow(deviatoric_eigen_values(2, 2), alpha_1)) / 3.0) + 
                      mu_2 * (std::pow(deviatoric_eigen_values(i, i), alpha_2) - (std::pow(deviatoric_eigen_values(0, 0), alpha_2) + std::pow(deviatoric_eigen_values(1, 1), alpha_2) + std::pow(deviatoric_eigen_values(2, 2), alpha_2)) / 3.0) +
                      mu_3 * (std::pow(deviatoric_eigen_values(i, i), alpha_3) - (std::pow(deviatoric_eigen_values(0, 0), alpha_3) + std::pow(deviatoric_eigen_values(1, 1), alpha_3) + std::pow(deviatoric_eigen_values(2, 2), alpha_3)) / 3.0);
        }

        // Compute the stress tensor
        BoundedMatrixType stress_tensor;
        noalias(stress_tensor) = beta(0) * M_container[0] + beta(1) * M_container[1] + beta(2) * M_container[2];


        // Stress vector computed
        noalias(r_stress_vector) = MathUtils<double>::StressTensorToVector(stress_tensor, VoigtSize);


        // We compute the tangent tensor by perturbation
        if (r_flags.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
            this->CalculateTangentTensor(rValues, ConstitutiveLaw::StressMeasure_PK2); // this modifies the ConstitutiveMatrix
        }
    }

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
    rFeatures.mOptions.Set(THREE_DIMENSIONAL_LAW);
    rFeatures.mOptions.Set(FINITE_STRAINS);
    rFeatures.mOptions.Set(ISOTROPIC);

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
    KRATOS_CHECK_VARIABLE_KEY(OGDEN_HYPERELASTIC_PARAMETERS);

    const Vector& r_ogden_props = rMaterialProperties[OGDEN_HYPERELASTIC_PARAMETERS];
    KRATOS_ERROR_IF(r_ogden_props.size() != VoigtSize) << "OGDEN_HYPERELASTIC_PARAMETERS must have 6 parameters" << std::endl;
    KRATOS_ERROR_IF(norm_2(r_ogden_props) <=  tolerance) << "OGDEN_HYPERELASTIC_PARAMETERS must have non-zero values" << std::endl;

    return 0;
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

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicOgden3D::CalculateTangentTensor(
    ConstitutiveLaw::Parameters& rValues,
    const ConstitutiveLaw::StressMeasure& rStressMeasure
    )
{
    // Calculates the Tangent Constitutive Tensor by perturbation
    TangentOperatorCalculatorUtility::CalculateTangentTensorFiniteDeformation(rValues, this, rStressMeasure);
}

/***********************************************************************************/
/***********************************************************************************/

} // Namespace Kratos
