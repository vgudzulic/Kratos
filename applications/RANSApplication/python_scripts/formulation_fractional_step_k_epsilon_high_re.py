from __future__ import absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics as Kratos
from KratosMultiphysics import IsDistributedRun
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable

# imports from fluid dynamics application
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import KratosMultiphysics.RANSApplication as KratosRANS

# imports from rans application
from KratosMultiphysics.RANSApplication.formulation import Formulation
from KratosMultiphysics.RANSApplication.strategy_factory import CreateStrategy
from KratosMultiphysics.RANSApplication.model_part_factory import CreateDuplicateModelPart
from KratosMultiphysics.RANSApplication.formulation_fractional_step import FractionalStepFormulation
from KratosMultiphysics.RANSApplication import RansCalculationUtilities

# case specific imports
if (IsDistributedRun() and CheckIfApplicationsAvailable("TrilinosApplication")):
    from KratosMultiphysics.RANSApplication.TrilinosExtension import MPICoupledStrategyItem as coupled_strategy_item
elif (not IsDistributedRun()):
    from KratosMultiphysics.RANSApplication import CoupledStrategyItem as coupled_strategy_item
else:
    raise Exception("Distributed run requires TrilinosApplication")


class FractionalStepKEpsilonHighReFormulation(Formulation):
    def __init__(self, model_part, settings):
        super(FractionalStepKEpsilonHighReFormulation, self).__init__(model_part, settings)
        self.defaults = Kratos.Parameters(r'''{
            "time_scheme_settings": {
                "scheme_type"         : "steady",
                "scheme_settings": {
                    "relative_velocity_tolerance": 1e-3,
                    "absolute_velocity_tolerance": 1e-5,
                    "relative_pressure_tolerance": 1e-3,
                    "absolute_pressure_tolerance": 1e-5,
                    "relative_turbulent_kinetic_energy_tolerance": 1e-3,
                    "absolute_turbulent_kinetic_energy_tolerance": 1e-5,
                    "relative_turbulent_energy_dissipation_rate_tolerance": 1e-3,
                    "absolute_turbulent_energy_dissipation_rate_tolerance": 1e-5,
                    "relative_turbulent_viscosity_tolerance": 1e-3,
                    "absolute_turbulent_viscosity_tolerance": 1e-5
                }
            },
            "constants":
            {
                "wall_smoothness_beta"                          : 5.2,
                "von_karman"                                    : 0.41,
                "c_mu"                                          : 0.09,
                "c1"                                            : 1.44,
                "c2"                                            : 1.92,
                "sigma_k"                                       : 1.0,
                "sigma_epsilon"                                 : 1.3,
                "stabilizing_upwind_operator_coefficient"       : 1.2,
                "stabilizing_positivity_preserving_coefficient" : 1.2
            },
            "consider_periodic_conditions": false,
            "velocity_pressure_coupling_iterations": 40,
            "turbulence_model_coupling_iterations": 20,
            "velocity_pressure_solver_settings": {},
            "turbulent_kinetic_energy_solver_settings": {
                "update_processes_list": [],
                "strategy_settings": {
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
                "update_processes_list": [],
                "strategy_settings": {
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

        self.settings.ValidateAndAssignDefaults(self.defaults)
        self.fractional_step_formulation = FractionalStepFormulation(model_part, settings["velocity_pressure_solver_settings"])
        self.SetIsPeriodic(self.settings["consider_periodic_conditions"].GetBool())
        self.fractional_step_formulation.SetIsPeriodic(self.IsPeriodic())

        Kratos.Logger.PrintInfo("FractionalStepKEpsilonHighReFormulation", "Construction of formulation finished.")

    def AddVariables(self):
        self.fractional_step_formulation.AddVariables()

        # RANS k - epsilon specific variables
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.RELAXED_ACCELERATION)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.KINEMATIC_VISCOSITY)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.TURBULENT_VISCOSITY)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(KratosRANS.TURBULENT_KINETIC_ENERGY)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(KratosRANS.TURBULENT_KINETIC_ENERGY_RATE)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE_2)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(KratosRANS.RANS_AUXILIARY_VARIABLE_1)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(KratosRANS.RANS_AUXILIARY_VARIABLE_2)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(KratosRANS.AFC_POSITIVE_ANTI_DIFFUSIVE_FLUX)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(KratosRANS.AFC_POSITIVE_ANTI_DIFFUSIVE_FLUX_LIMIT)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(KratosRANS.AFC_NEGATIVE_ANTI_DIFFUSIVE_FLUX)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(KratosRANS.AFC_NEGATIVE_ANTI_DIFFUSIVE_FLUX_LIMIT)

        Kratos.Logger.PrintInfo("FractionalStepKEpsilonHighReFormulation", "Added solution step variables.")

    def AddDofs(self):
        self.fractional_step_formulation.AddDofs()

        Kratos.VariableUtils().AddDof(KratosRANS.TURBULENT_KINETIC_ENERGY, self.GetBaseModelPart())
        Kratos.VariableUtils().AddDof(KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE, self.GetBaseModelPart())

        Kratos.Logger.PrintInfo("FractionalStepKEpsilonHighReFormulation", "Added solution dofs.")

    def PrepareModelPart(self):
        self.domain_size = self.GetBaseModelPart().ProcessInfo[Kratos.DOMAIN_SIZE]
        condition_suffix = str(self.domain_size) + "D" + str(self.domain_size) + "N"
        element_suffix = str(self.domain_size) + "D" + str(self.domain_size + 1) + "N"

        self.fractional_step_formulation.PrepareModelPart()

        # preparing tke model part
        element_name = "RansEvmKEpsilonK" + element_suffix
        condition_name = "Condition" + condition_suffix
        modelpart_name = self.GetBaseModelPart().Name + "_" + element_name + "_" + condition_name
        self.tke_model_part = CreateDuplicateModelPart(self.GetBaseModelPart(), modelpart_name,
                                                            element_name, condition_name, "")

        # preparing epsilon model part
        element_name = "RansEvmKEpsilonEpsilon" + element_suffix
        condition_name = "RansEvmKEpsilonEpsilonWall" + condition_suffix
        modelpart_name = self.GetBaseModelPart().Name + "_" + element_name + "_" + condition_name
        self.epsilon_model_part = CreateDuplicateModelPart(self.GetBaseModelPart(), modelpart_name,
                                                            element_name, condition_name, "")

        Kratos.Logger.PrintInfo("FractionalStepKEpsilonHighReFormulation", "Solver model parts prepared.")

    def Initialize(self):
        self.__InitializeConstants()
        self.fractional_step_formulation.Initialize()

        time_scheme_settings = self.settings["time_scheme_settings"]

        self.tke_strategy = CreateStrategy(
                                    self.settings["turbulent_kinetic_energy_solver_settings"]["strategy_settings"],
                                    time_scheme_settings,
                                    self.tke_model_part,
                                    self.GetBaseModelPart(),
                                    self.GetCommunicator(),
                                    KratosRANS.TURBULENT_KINETIC_ENERGY,
                                    KratosRANS.TURBULENT_KINETIC_ENERGY_RATE,
                                    KratosRANS.RANS_AUXILIARY_VARIABLE_1)

        self.epsilon_strategy = CreateStrategy(
                                    self.settings["turbulent_energy_dissipation_rate_solver_settings"]["strategy_settings"],
                                    time_scheme_settings,
                                    self.epsilon_model_part,
                                    self.GetBaseModelPart(),
                                    self.GetCommunicator(),
                                    KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE,
                                    KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE_2,
                                    KratosRANS.RANS_AUXILIARY_VARIABLE_2)

        time_scheme_type = time_scheme_settings["scheme_type"].GetString()

        if (time_scheme_type == "steady"):
            self.is_steady = True
            tke_strategy_item = coupled_strategy_item(
                                    self.tke_strategy,
                                    "TurbulentKineticEnergy",
                                    self.settings["turbulent_kinetic_energy_solver_settings"]["strategy_settings"]["echo_level"].GetInt())
            epsilon_strategy_item = coupled_strategy_item(
                                    self.epsilon_strategy,
                                    "TurbulentEnergyDissipationRate",
                                    self.settings["turbulent_energy_dissipation_rate_solver_settings"]["strategy_settings"]["echo_level"].GetInt()
                                    )
        elif (time_scheme_type == "transient"):
            self.is_steady = False
            tke_strategy_item = coupled_strategy_item(
                                    self.tke_strategy,
                                    "TurbulentKineticEnergy",
                                    self.settings["turbulent_kinetic_energy_solver_settings"]["echo_level"].GetInt())
            epsilon_strategy_item = coupled_strategy_item(
                                    self.tke_strategy,
                                    "TurbulentKineticEnergy",
                                    self.settings["turbulent_energy_dissipation_rate_solver_settings"]["echo_level"].GetInt())

        velocity_pressure_coupling_iterations = self.settings["velocity_pressure_coupling_iterations"].GetInt()
        turbulence_model_coupling_iterations = self.settings["turbulence_model_coupling_iterations"].GetInt()

        self.velocity_pressure_solvability_pattern = [1] * velocity_pressure_coupling_iterations
        self.velocity_pressure_solvability_pattern.extend([0] * turbulence_model_coupling_iterations)

        self.turbulence_model_solvability_pattern = [0] * velocity_pressure_coupling_iterations
        self.turbulence_model_solvability_pattern.extend([1] * turbulence_model_coupling_iterations)

        self.strategy_items = self.fractional_step_formulation.GetCoupledStrategyItems()
        # self.strategy_items[-1][1].SetStrategySolvabilityPattern(self.velocity_pressure_solvability_pattern)
        self.strategy_items.append(["turbulent_kinetic_energy_solver_settings", tke_strategy_item])
        # self.strategy_items[-1][1].SetStrategySolvabilityPattern(self.turbulence_model_solvability_pattern)
        self.strategy_items.append(["turbulent_energy_dissipation_rate_solver_settings", epsilon_strategy_item])
        # self.strategy_items[-1][1].SetStrategySolvabilityPattern(self.turbulence_model_solvability_pattern)

        Kratos.Logger.PrintInfo("FractionalStepKEpsilonHighReFormulation", "Solver initialization finished.")

    def __InitializeConstants(self):
        # update model constants
        constants = self.settings["constants"]
        constants.ValidateAndAssignDefaults(self.defaults["constants"])
        self.GetBaseModelPart().ProcessInfo[KratosRANS.WALL_SMOOTHNESS_BETA] = constants["wall_smoothness_beta"].GetDouble()
        self.GetBaseModelPart().ProcessInfo[KratosRANS.WALL_VON_KARMAN] = constants["von_karman"].GetDouble()
        self.GetBaseModelPart().ProcessInfo[KratosRANS.TURBULENCE_RANS_C_MU] = constants["c_mu"].GetDouble()
        self.GetBaseModelPart().ProcessInfo[KratosRANS.TURBULENCE_RANS_C1] = constants["c1"].GetDouble()
        self.GetBaseModelPart().ProcessInfo[KratosRANS.TURBULENCE_RANS_C2] = constants["c2"].GetDouble()
        self.GetBaseModelPart().ProcessInfo[KratosRANS.TURBULENT_KINETIC_ENERGY_SIGMA] = constants["sigma_k"].GetDouble()
        self.GetBaseModelPart().ProcessInfo[KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA] = constants["sigma_epsilon"].GetDouble()
        self.GetBaseModelPart().ProcessInfo[KratosRANS.RANS_STABILIZATION_DISCRETE_UPWIND_OPERATOR_COEFFICIENT] = constants["stabilizing_upwind_operator_coefficient"].GetDouble()
        self.GetBaseModelPart().ProcessInfo[KratosRANS.RANS_STABILIZATION_DIAGONAL_POSITIVITY_PRESERVING_COEFFICIENT] = constants["stabilizing_positivity_preserving_coefficient"].GetDouble()
        self.GetBaseModelPart().ProcessInfo[KratosRANS.RANS_Y_PLUS_LIMIT] = RansCalculationUtilities.CalculateLogarithmicYPlusLimit(
                                                                                self.GetBaseModelPart().ProcessInfo[KratosRANS.WALL_VON_KARMAN],
                                                                                self.GetBaseModelPart().ProcessInfo[KratosRANS.WALL_SMOOTHNESS_BETA]
                                                                                )
    def GetCoupledStrategyItems(self):
        return self.strategy_items

    def Finalize(self):
        self.fractional_step_formulation.Finalize()
        self.tke_strategy.Clear()
        self.epsilon_strategy.Clear()

    def GetMinimumBufferSize(self):
        return 3

    def IsSolvingForSteadyState(self):
        return self.is_steady

    def SetCommunicator(self, communicator):
        super(FractionalStepKEpsilonHighReFormulation, self).SetCommunicator(communicator)
        self.fractional_step_formulation.SetCommunicator(communicator)