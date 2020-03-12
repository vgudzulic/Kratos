from abc import ABC, abstractmethod

# imports from kratos core
import KratosMultiphysics as Kratos
from KratosMultiphysics.process_factory import KratosProcessFactory

# imports from fluid dynamics application
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

# imports from rans application
import KratosMultiphysics.RANSApplication as KratosRANS
from KratosMultiphysics.RANSApplication.model_part_factory import CreateDuplicateModelPart
from KratosMultiphysics.RANSApplication.strategy_factory import CreateStrategy
from KratosMultiphysics.RANSApplication import RansCalculationUtilities

def CreateFormulation(model_part, formulation_name, settings, scheme_settings):
    formulations_list = [
        ["segregated_vms_k_epsilon_high_re", SegregatedVMSKEpsilonHighRe]
    ]

    formulation_names_list = [formulations_list[i][0] for i in range(len(formulations_list))]
    formulation_list = [formulations_list[i][1] for i in range(len(formulations_list))]

    if (formulation_name not in formulation_names_list):
        msg = "Unknown formulation type=\"" + formulation_name + "\". \nFollowing formulations are allowed:\n    "
        msg += "\n    ".join(sorted(formulation_names_list))
        raise Exception(msg + "\n")

    current_formulation = formulation_list[formulation_names_list.index(formulation_name)](model_part, settings, scheme_settings)

    Kratos.Logger.PrintInfo("RANSFormulationFactory", "Created " + formulation_name + " formulation.")

    return current_formulation

class Formulation(ABC):
    def __init__(self, base_model_part, coupling_settings, scheme_settings):
        self.settings = coupling_settings
        self.scheme_settings = scheme_settings
        self.base_model_part = base_model_part
        self.strategy_list = []
        self.communicator = None
        self.solver_names_list = []
        self.variable_list = []
        self.model_part_list = []

        defaults = Kratos.Parameters(r'''{

        }''')

    @abstractmethod
    def AddVariables(self):
        pass

    @abstractmethod
    def AddDofs(self):
        pass

    @abstractmethod
    def PrepareModelPart(self):
        pass

    def Initialize(self):
        if self.scheme_settings["scheme_type"].GetString() != "steady":
            raise Exception("Currently only steady schemes are supported")

        default_strategy_settings = Kratos.Parameters(r'''{
            "relative_tolerance": 1e-3,
            "absolute_tolerance": 1e-5,
            "update_processes": [],
            "strategy_settings": {}
        }''')

        for i, solver_name in enumerate(self.solver_names_list):
            settings = self.settings[solver_name]
            settings.ValidateAndAssignDefaults(default_strategy_settings)
            solver_strategy = {}
            solver_strategy["strategy"] = CreateStrategy(settings["strategy_settings"],
                                                         self.scheme_settings,
                                                         self.model_part_list[i],
                                                         self.base_model_part,
                                                         self.communicator,
                                                         self.variable_list[i][0],
                                                         self.variable_list[i][1],
                                                         self.variable_list[i][2])

            solver_strategy["relative_tolerance"] = settings["relative_tolerance"].GetDouble()
            solver_strategy["absolute_tolerance"] = settings["absolute_tolerance"].GetDouble()
            solver_strategy["variable_name"] = self.variable_list[i][0].Name()

            factory = KratosProcessFactory(self.base_model_part.GetModel())
            solver_strategy["update_processes_list"] = factory.ConstructListOfProcesses(settings["update_processes"])

            self.strategy_list.append(solver_strategy)

    def SetCommunicator(self, communicator):
        self.communicator = communicator

    def GetBaseModelPart(self):
        return self.base_model_part

    def GetStrategyList(self):
        return self.strategy_list

    def IsPeriodic(self):
        is_periodic = False
        for solver_name in self.solver_names_list:
            if not is_periodic: is_periodic = self.settings[solver_name]["strategy_settings"]["is_periodic"].GetBool()


class SegregatedVMSKEpsilonHighRe(Formulation):
    def __init__(self, model, settings, scheme_settings):
        super(SegregatedVMSKEpsilonHighRe, self).__init__(model, settings, scheme_settings)
        defaults = Kratos.Parameters(r'''{
            "constants":
            {
                "wall_smoothness_beta"    : 5.2,
                "von_karman"              : 0.41,
                "c_mu"                    : 0.09,
                "c1"                      : 1.44,
                "c2"                      : 1.92,
                "sigma_k"                 : 1.0,
                "sigma_epsilon"           : 1.3
            },
            "velocity_solver_settings": {
                "relative_tolerance": 1e-3,
                "absolute_tolerance": 1e-5,
                "update_processes": [],
                "strategy_settings": {
                    "is_periodic"           : false,
                    "relative_tolerance"    : 1e-3,
                    "absolute_tolerance"    : 1e-5,
                    "max_iterations"        : 200,
                    "relaxation_factor"     : 0.5,
                    "echo_level"            : 2,
                    "linear_solver_settings": {
                        "solver_type"  : "amgcl"
                    },
                    "reform_dofs_at_each_step": true,
                    "move_mesh_strategy": 0,
                    "move_mesh_flag": false,
                    "compute_reactions": false
                }
            },
            "pressure_solver_settings": {
                "relative_tolerance": 1e-3,
                "absolute_tolerance": 1e-5,
                "update_processes": [],
                "strategy_settings": {
                    "is_periodic"           : false,
                    "relative_tolerance"    : 1e-3,
                    "absolute_tolerance"    : 1e-5,
                    "max_iterations"        : 200,
                    "relaxation_factor"     : 0.5,
                    "echo_level"            : 2,
                    "linear_solver_settings": {
                        "solver_type"  : "amgcl"
                    },
                    "reform_dofs_at_each_step": true,
                    "move_mesh_strategy": 0,
                    "move_mesh_flag": false,
                    "compute_reactions": false
                }
            },
            "turbulent_kinetic_energy_solver_settings": {
                "relative_tolerance": 1e-3,
                "absolute_tolerance": 1e-5,
                "update_processes": [],
                "strategy_settings": {
                    "is_periodic"           : false,
                    "relative_tolerance"    : 1e-3,
                    "absolute_tolerance"    : 1e-5,
                    "max_iterations"        : 200,
                    "relaxation_factor"     : 0.5,
                    "echo_level"            : 2,
                    "linear_solver_settings": {
                        "solver_type"  : "amgcl"
                    },
                    "reform_dofs_at_each_step": true,
                    "move_mesh_strategy": 0,
                    "move_mesh_flag": false,
                    "compute_reactions": false
                }
            },
            "turbulent_energy_dissipation_rate_solver_settings": {
                "relative_tolerance": 1e-3,
                "absolute_tolerance": 1e-5,
                "update_processes": [],
                "strategy_settings": {
                    "is_periodic"           : false,
                    "relative_tolerance"    : 1e-3,
                    "absolute_tolerance"    : 1e-5,
                    "max_iterations"        : 200,
                    "relaxation_factor"     : 0.5,
                    "echo_level"            : 2,
                    "linear_solver_settings": {
                        "solver_type"  : "amgcl"
                    },
                    "reform_dofs_at_each_step": true,
                    "move_mesh_strategy": 0,
                    "move_mesh_flag": false,
                    "compute_reactions": false
                }
            }
        }''')

        self.settings.RecursivelyAddMissingParameters(defaults)

        self.solver_names_list = [
            "velocity_solver_settings",
            "pressure_solver_settings",
            "turbulent_kinetic_energy_solver_settings",
            "turbulent_energy_dissipation_rate_solver_settings"
            ]

        self.variable_list = [
            [Kratos.VELOCITY, Kratos.ACCELERATION, Kratos.RELAXED_ACCELERATION],
            [Kratos.PRESSURE, None, None],
            [KratosRANS.TURBULENT_KINETIC_ENERGY, KratosRANS.TURBULENT_KINETIC_ENERGY_RATE, KratosRANS.RANS_AUXILIARY_VARIABLE_1],
            [KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE, KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE_2, KratosRANS.RANS_AUXILIARY_VARIABLE_2]
        ]

    def AddVariables(self):
        self.base_model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        self.base_model_part.AddNodalSolutionStepVariable(Kratos.ACCELERATION)
        self.base_model_part.AddNodalSolutionStepVariable(Kratos.RELAXED_ACCELERATION)
        self.base_model_part.AddNodalSolutionStepVariable(Kratos.MESH_VELOCITY)
        self.base_model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        self.base_model_part.AddNodalSolutionStepVariable(Kratos.EXTERNAL_PRESSURE)
        self.base_model_part.AddNodalSolutionStepVariable(Kratos.VISCOSITY)
        self.base_model_part.AddNodalSolutionStepVariable(Kratos.DENSITY)
        self.base_model_part.AddNodalSolutionStepVariable(Kratos.BODY_FORCE)
        self.base_model_part.AddNodalSolutionStepVariable(Kratos.REACTION)
        self.base_model_part.AddNodalSolutionStepVariable(Kratos.REACTION_WATER_PRESSURE)
        self.base_model_part.AddNodalSolutionStepVariable(Kratos.NORMAL)
        self.base_model_part.AddNodalSolutionStepVariable(Kratos.KINEMATIC_VISCOSITY)
        self.base_model_part.AddNodalSolutionStepVariable(Kratos.TURBULENT_VISCOSITY)
        self.base_model_part.AddNodalSolutionStepVariable(KratosRANS.TURBULENT_KINETIC_ENERGY)
        self.base_model_part.AddNodalSolutionStepVariable(KratosRANS.TURBULENT_KINETIC_ENERGY_RATE)
        self.base_model_part.AddNodalSolutionStepVariable(KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE)
        self.base_model_part.AddNodalSolutionStepVariable(KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE_2)
        self.base_model_part.AddNodalSolutionStepVariable(KratosRANS.RANS_AUXILIARY_VARIABLE_1)
        self.base_model_part.AddNodalSolutionStepVariable(KratosRANS.RANS_AUXILIARY_VARIABLE_2)

        if self.IsPeriodic():
            self.base_model_part.AddNodalSolutionStepVariable(KratosCFD.PATCH_INDEX)

    def AddDofs(self):
        Kratos.VariableUtils().AddDof(Kratos.VELOCITY_X, Kratos.REACTION_X,self.base_model_part)
        Kratos.VariableUtils().AddDof(Kratos.VELOCITY_Y, Kratos.REACTION_Y,self.base_model_part)
        Kratos.VariableUtils().AddDof(Kratos.VELOCITY_Z, Kratos.REACTION_Z,self.base_model_part)
        Kratos.VariableUtils().AddDof(Kratos.PRESSURE, Kratos.REACTION_WATER_PRESSURE,self.base_model_part)
        Kratos.VariableUtils().AddDof(KratosRANS.TURBULENT_KINETIC_ENERGY, self.base_model_part)
        Kratos.VariableUtils().AddDof(KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE, self.base_model_part)

    def PrepareModelPart(self):
        domain_size = self.base_model_part.ProcessInfo[Kratos.DOMAIN_SIZE]
        condition_suffix = str(domain_size) + "D" + str(domain_size) + "N"
        element_suffix = str(domain_size) + "D" + str(domain_size + 1) + "N"

        # preparing velocity model part
        element_name = "SegregatedVMSVelocity" + element_suffix
        condition_name = "SegregatedVMSWallVelocity" + condition_suffix
        modelpart_name = self.base_model_part.Name + "_" + element_name + "_" + condition_name
        self.model_part_list.append(CreateDuplicateModelPart(self.base_model_part, modelpart_name,
                                                            element_name, condition_name, ""))

        # preparing pressure model part
        element_name = "SegregatedVMSPressure" + element_suffix
        condition_name = "Condition" + condition_suffix
        modelpart_name = self.base_model_part.Name + "_" + element_name + "_" + condition_name
        self.model_part_list.append(CreateDuplicateModelPart(self.base_model_part, modelpart_name,
                                                            element_name, condition_name, ""))

        # preparing tke model part
        element_name = "RansEvmKEpsilonK" + element_suffix
        condition_name = "Condition" + condition_suffix
        modelpart_name = self.base_model_part.Name + "_" + element_name + "_" + condition_name
        self.model_part_list.append(CreateDuplicateModelPart(self.base_model_part, modelpart_name,
                                                            element_name, condition_name, ""))

        # preparing epsilon model part
        element_name = "RansEvmKEpsilonEpsilon" + element_suffix
        condition_name = "RansEvmKEpsilonEpsilonWall" + condition_suffix
        modelpart_name = self.base_model_part.Name + "_" + element_name + "_" + condition_name
        self.model_part_list.append(CreateDuplicateModelPart(self.base_model_part, modelpart_name,
                                                            element_name, condition_name, ""))

    def Initialize(self):
        # update model constants
        constants = self.settings["constants"]
        self.base_model_part.ProcessInfo[KratosRANS.WALL_SMOOTHNESS_BETA] = constants["wall_smoothness_beta"].GetDouble()
        self.base_model_part.ProcessInfo[KratosRANS.WALL_VON_KARMAN] = constants["von_karman"].GetDouble()
        self.base_model_part.ProcessInfo[KratosRANS.TURBULENCE_RANS_C_MU] = constants["c_mu"].GetDouble()
        self.base_model_part.ProcessInfo[KratosRANS.TURBULENCE_RANS_C1] = constants["c1"].GetDouble()
        self.base_model_part.ProcessInfo[KratosRANS.TURBULENCE_RANS_C2] = constants["c2"].GetDouble()
        self.base_model_part.ProcessInfo[KratosRANS.TURBULENT_KINETIC_ENERGY_SIGMA] = constants["sigma_k"].GetDouble()
        self.base_model_part.ProcessInfo[KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA] = constants["sigma_epsilon"].GetDouble()
        self.base_model_part.ProcessInfo[KratosRANS.RANS_STABILIZATION_DISCRETE_UPWIND_OPERATOR_COEFFICIENT] = 1.2
        self.base_model_part.ProcessInfo[KratosRANS.RANS_STABILIZATION_DIAGONAL_POSITIVITY_PRESERVING_COEFFICIENT] = 1.2
        self.base_model_part.ProcessInfo[KratosRANS.RANS_Y_PLUS_LIMIT] = RansCalculationUtilities.CalculateLogarithmicYPlusLimit(
                                                                                self.base_model_part.ProcessInfo[KratosRANS.WALL_VON_KARMAN],
                                                                                self.base_model_part.ProcessInfo[KratosRANS.WALL_SMOOTHNESS_BETA]
                                                                                )
        super(SegregatedVMSKEpsilonHighRe, self).Initialize()

