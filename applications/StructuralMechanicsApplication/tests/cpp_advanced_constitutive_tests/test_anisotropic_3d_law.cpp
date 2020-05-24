// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                     license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
// #include "includes/gid_io.h"
#include "utilities/read_materials_utility.h"
#include "custom_advanced_constitutive/generic_anisotropic_3d_law.h"
#include "includes/mat_variables.h"

namespace Kratos
{
namespace Testing
{
/// Nodetype definition
typedef Node<3> NodeType;

// void GiDIODebugAnisotropic(ModelPart& ThisModelPart)
// {
//     GidIO<> gid_io("TEST_ANISOTROPIC", GiD_PostBinary, SingleFile, WriteUndeformed,  WriteElementsOnly);
//     const int nl_iter = ThisModelPart.GetProcessInfo()[NL_ITERATION_NUMBER];
//     const double label = static_cast<double>(nl_iter);
//
//     gid_io.InitializeMesh(label);
//     gid_io.WriteMesh(ThisModelPart.GetMesh());
//     gid_io.FinalizeMesh();
//     gid_io.InitializeResults(label, ThisModelPart.GetMesh());
//     gid_io.WriteNodalResults(DISPLACEMENT, ThisModelPart.Nodes(), label, 0);
//     gid_io.PrintOnGaussPoints(GREEN_LAGRANGE_STRAIN_VECTOR, ThisModelPart, label);
//     gid_io.PrintOnGaussPoints(PK2_STRESS_VECTOR, ThisModelPart, label);
// }

void Create3DGeometryHexahedraAnisotropicLaw(ModelPart& rThisModelPart, const std::string ElementName = "SmallDisplacementElement3D8N")
{
    rThisModelPart.AddNodalSolutionStepVariable(DISPLACEMENT);

    auto& r_process_info = rThisModelPart.GetProcessInfo();
    r_process_info[STEP] = 1;
    r_process_info[NL_ITERATION_NUMBER] = 1;

    Parameters parameters = Parameters(R"(
    {
    "properties" : [{
        "model_part_name" : "Main",
        "properties_id"   : 1,
        "Material"        : {
            "constitutive_law" : {
                "name" : "GenericAnisotropic3DLaw"
            },
            "Variables"        : {
                "ORTHOTROPIC_ELASTIC_CONSTANTS"     : [40e9, 10e9, 10e9, 0.2, 0.2, 0.2],
                "EULER_ANGLES"                      : [45,0,0],
                "ISOTROPIC_ANISOTROPIC_YIELD_RATIO" : [2.0,0.667,0.667,1.0,1.0,1.0]
            },
            "Tables"           : {}
        },
        "sub_properties" : [{
            "properties_id"   : 10,
            "Material"        : {
                "constitutive_law" : {
                    "name" : "LinearElastic3DLaw"
                },
                "Variables"        : {
                    "DENSITY"         : 2400.0,
                    "YOUNG_MODULUS"   : 40E9,
                    "POISSON_RATIO"   : 0.3
                },
                "Tables"           : {}
            }
        }]
    }]
    })" );

    // Read properties
    auto read_util = ReadMaterialsUtility(rThisModelPart.GetModel());
    read_util.ReadMaterials(parameters);

    // Create nodes and elements
    Properties::Pointer p_elem_prop = rThisModelPart.pGetProperties(1);

    // First we create the nodes
    rThisModelPart.CreateNewNode(1 , 0.0 , 1.0 , 1.0);
    rThisModelPart.CreateNewNode(2 , 0.0 , 1.0 , 0.0);
    rThisModelPart.CreateNewNode(3 , 0.0 , 0.0 , 1.0);
    rThisModelPart.CreateNewNode(4 , 1.0 , 1.0 , 1.0);
    rThisModelPart.CreateNewNode(5 , 0.0 , 0.0 , 0.0);
    rThisModelPart.CreateNewNode(6 , 1.0 , 1.0 , 0.0);
    rThisModelPart.CreateNewNode(7 , 1.0 , 0.0 , 1.0);
    rThisModelPart.CreateNewNode(8 , 1.0 , 0.0 , 0.0);

    // Now we create the elements
    rThisModelPart.CreateNewElement(ElementName, 1, {{5,8,6,2,3,7,4,1}}, p_elem_prop);

    const auto& const_process_info = rThisModelPart.GetProcessInfo();

    // Initialize elements
    for (auto& r_elem : rThisModelPart.Elements()) {
        r_elem.Initialize(const_process_info);
        r_elem.InitializeSolutionStep(const_process_info);
        r_elem.InitializeNonLinearIteration(const_process_info);
    }
}

/**
* Check the correct work of anisotropic 3D law
*/
KRATOS_TEST_CASE_IN_SUITE(AnisotropicConstitutiveLawHexahedron, KratosStructuralMechanicsFastSuite)
{
    Model this_model;
    ModelPart& r_model_part = this_model.CreateModelPart("Main", 2);
    Create3DGeometryHexahedraAnisotropicLaw(r_model_part);

    const array_1d<double, 3> zero = ZeroVector(3);
    array_1d<double, 3> delta = ZeroVector(3);
    delta[1] = 1.0e-3;
    for (auto& node : r_model_part.Nodes()) {
        if (node.X() < 1.0e-3) {
            node.Fix(DISPLACEMENT_X);
            node.Fix(DISPLACEMENT_Y);
            node.Fix(DISPLACEMENT_Z);
            node.FastGetSolutionStepValue(DISPLACEMENT) = zero;
        } else {
            node.FastGetSolutionStepValue(DISPLACEMENT) = delta;
            node.Coordinates() += delta;
        }
    }

//     // DEBUG
//     GiDIODebugAnisotropic(model_part);

    /// Tolerance
    const double tolerance = 1.0e-4;

    ProcessInfo& process_info = r_model_part.GetProcessInfo();
    for (auto& elem : r_model_part.Elements()) {
        std::vector<Vector> solution;
        elem.CalculateOnIntegrationPoints(PK2_STRESS_VECTOR, solution, process_info);

        Vector reference = ZeroVector(6);
        reference[0] = 7.61218e+06;
        reference[1] = 7.61218e+06;
        reference[2] = 160256;
        reference[3] = 1.16186e+07;

        for (auto& sol : solution) {
            KRATOS_CHECK_VECTOR_RELATIVE_NEAR(sol, reference, tolerance);
        }
    }
}

} // namespace Testing
} // namespace Kratos
