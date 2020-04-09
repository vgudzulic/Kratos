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

// System includes

// External includes

// Project includes
#include "includes/cfd_variables.h"
#include "includes/define.h"
#include "utilities/variable_utils.h"

// Application includes
#include "custom_utilities/rans_calculation_utilities.h"
#include "custom_utilities/rans_check_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "rans_wall_function_update_process.h"

namespace Kratos
{
RansWallFunctionUpdateProcess::RansWallFunctionUpdateProcess(Model& rModel, Parameters rParameters)
    : mrModel(rModel), mrParameters(rParameters)
{
    KRATOS_TRY

    Parameters default_parameters = Parameters(R"(
        {
            "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "echo_level"      : 0,
            "von_karman"      : 0.41,
            "beta"            : 5.2,
            "c_mu"            : 0.09
        })");

    mrParameters.ValidateAndAssignDefaults(default_parameters);

    mEchoLevel = mrParameters["echo_level"].GetInt();
    mModelPartName = mrParameters["model_part_name"].GetString();
    mVonKarman = mrParameters["von_karman"].GetDouble();
    mBeta = mrParameters["beta"].GetDouble();
    mCmu = mrParameters["c_mu"].GetDouble();
    mLimitYPlus =
        RansCalculationUtilities::CalculateLogarithmicYPlusLimit(mVonKarman, mBeta);

    KRATOS_CATCH("");
}

int RansWallFunctionUpdateProcess::Check()
{
    KRATOS_TRY

    RansCheckUtilities::CheckIfModelPartExists(mrModel, mModelPartName);

    const ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);

    RansCheckUtilities::CheckIfVariableExistsInModelPart(r_model_part, KINEMATIC_VISCOSITY);
    RansCheckUtilities::CheckIfVariableExistsInModelPart(r_model_part, TURBULENT_VISCOSITY);

    return 0;

    KRATOS_CATCH("");
}

void RansWallFunctionUpdateProcess::ExecuteInitialize()
{
}

void RansWallFunctionUpdateProcess::Execute()
{
    KRATOS_TRY

    ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);

    ModelPart::ConditionsContainerType& r_conditions = r_model_part.Conditions();
    const int number_of_conditions = r_conditions.size();

    const double c_mu_25 = std::pow(mCmu, 0.25);

#pragma omp parallel for
    for (int i_cond = 0; i_cond < number_of_conditions; ++i_cond)
    {
        ModelPart::ConditionType& r_condition = *(r_conditions.begin() + i_cond);

        const array_1d<double, 3> wall_cell_center_velocity =
            RansCalculationUtilities::CalculateWallVelocity(r_condition);
        const double wall_cell_center_velocity_magnitude = norm_2(wall_cell_center_velocity);
        const array_1d<double, 3>& r_normal = r_condition.GetValue(NORMAL);

        const double wall_height =
            RansCalculationUtilities::CalculateWallHeight(r_condition, r_normal);

        const double nu = RansCalculationUtilities::EvaluateInParentCenter(
            KINEMATIC_VISCOSITY, r_condition);

        double y_plus{0.0}, u_tau{0.0};
        RansCalculationUtilities::CalculateYPlusAndUtau(
            y_plus, u_tau, wall_cell_center_velocity_magnitude, wall_height, nu,
            mVonKarman, mBeta);

        // if (y_plus < mLimitYPlus)
        // {
        //     const ConditionType::GeometryType& r_geometry = r_condition.GetGeometry();
        //     Matrix shape_functions;
        //     Vector gauss_weights;
        //     RansCalculationUtilities::CalculateConditionGeometryData(
        //         r_geometry, GeometryData::GI_GAUSS_1, gauss_weights, shape_functions);

        //     const Vector& r_gauss_shape_functions = row(shape_functions, 0);

        //     const array_1d<double, 3>& r_wall_velocity =
        //         RansCalculationUtilities::EvaluateInPoint(
        //             r_geometry, VELOCITY, r_gauss_shape_functions);
        //     const double tke = RansCalculationUtilities::EvaluateInPoint(
        //         r_geometry, TURBULENT_KINETIC_ENERGY, r_gauss_shape_functions);
        //     const double wall_velocity_magnitude = norm_2(r_wall_velocity);
        //     u_tau = std::max(c_mu_25 * std::sqrt(std::max(tke, 0.0)),
        //                      wall_velocity_magnitude / mLimitYPlus);
        //     y_plus = u_tau * wall_height / nu;
        // }

        r_condition.SetValue(RANS_Y_PLUS, y_plus);
        r_condition.SetValue(FRICTION_VELOCITY, u_tau);
    }

    KRATOS_CATCH("");
}

std::string RansWallFunctionUpdateProcess::Info() const
{
    return std::string("RansWallFunctionUpdateProcess");
}

void RansWallFunctionUpdateProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

void RansWallFunctionUpdateProcess::PrintData(std::ostream& rOStream) const
{
}

} // namespace Kratos.
