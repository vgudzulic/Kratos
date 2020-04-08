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
#include "input_output/vtk_output.h"

// Application includes
#include "rans_application_variables.h"

// Include base h
#include "rans_vtk_output_process.h"

namespace Kratos
{
RansVtkOutputProcess::RansVtkOutputProcess(Model& rModel, Parameters rParameters)
    : mrModel(rModel), mrParameters(rParameters)
{
    KRATOS_TRY

    Parameters default_parameters = Parameters(R"(
        {
            "model_part_name"          : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "output_file_name_pattern" : "NL_<model_part_name>_R<rank>_S<step>_CL<coupling_iteration>_NL<nl_iteration>"
        })");

    mrParameters.AddMissingParameters(default_parameters);
    mOutputFilenamePattern = mrParameters["output_file_name_pattern"].GetString();
    mrParameters.RemoveValue("output_file_name_pattern");

    ModelPart& r_model_part =
        mrModel.GetModelPart(mrParameters["model_part_name"].GetString());

    mpVtkOutput = Kratos::make_unique<VtkOutput>(r_model_part, mrParameters);

    KRATOS_CATCH("");
}

void RansVtkOutputProcess::Execute()
{
    std::string current_filename = mOutputFilenamePattern;

    const std::string model_part_name = mrParameters["model_part_name"].GetString();
    ModelPart& r_model_part = mrModel.GetModelPart(model_part_name);

    const ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
    const int step = r_process_info[STEP];
    const int nl_iteration = r_process_info[NL_ITERATION_NUMBER];
    const int coupling_iteration = r_process_info[COUPLING_ITERATION];
    const int rank = r_model_part.GetCommunicator().MyPID();

    ReplaceAll(current_filename, "<model_part_name>", model_part_name);
    ReplaceAll(current_filename, "<rank>", std::to_string(rank));
    ReplaceAll(current_filename, "<step>", std::to_string(step));
    ReplaceAll(current_filename, "<nl_iteration>", std::to_string(nl_iteration));
    ReplaceAll(current_filename, "<coupling_iteration>", std::to_string(coupling_iteration));

    mpVtkOutput->PrintOutput(current_filename);
}

void RansVtkOutputProcess::ReplaceAll(std::string& rOutput,
                                      const std::string& rFindString,
                                      const std::string& rReplaceString)
{
    // Get the first occurrence
    size_t pos = rOutput.find(rFindString);

    // Repeat till end is reached
    while (pos != std::string::npos)
    {
        // Replace this occurrence of Sub String
        rOutput.replace(pos, rFindString.size(), rReplaceString);
        // Get the next occurrence from the current position
        pos = rOutput.find(rFindString, pos + rReplaceString.size());
    }
}

std::string RansVtkOutputProcess::Info() const
{
    return std::string("RansVtkOutputProcess");
}

void RansVtkOutputProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

void RansVtkOutputProcess::PrintData(std::ostream& rOStream) const
{
}
} // namespace Kratos.