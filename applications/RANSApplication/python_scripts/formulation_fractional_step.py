from __future__ import absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics as Kratos
from KratosMultiphysics import IsDistributedRun
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

# application imports
from KratosMultiphysics.RANSApplication.formulation import Formulation
from KratosMultiphysics.RANSApplication.strategy_factory import CreateLinearSolver
from KratosMultiphysics.RANSApplication.model_part_factory import CreateDuplicateModelPart

# case specific imports
if (IsDistributedRun() and CheckIfApplicationsAvailable("TrilinosApplication")):
    from KratosMultiphysics.FluidDynamicsApplication.TrilinosExtension import TrilinosFractionalStepSettingsPeriodic as periodic_solver_settings
    from KratosMultiphysics.FluidDynamicsApplication.TrilinosExtension import TrilinosFractionalStepSettings as solver_settings
    from KratosMultiphysics.FluidDynamicsApplication.TrilinosExtension import TrilinosFSStrategy as solver_strategy
    from KratosMultiphysics.FluidDynamicsApplication.TrilinosExtension import TrilinosStrategyLabel as strategy_label
    from KratosMultiphysics.RANSApplication.TrilinosExtension import MPICoupledStrategyItem as coupled_strategy_item
elif (not IsDistributedRun()):
    from KratosMultiphysics.FluidDynamicsApplication import FractionalStepSettingsPeriodic as periodic_solver_settings
    from KratosMultiphysics.FluidDynamicsApplication import FractionalStepSettings as solver_settings
    from KratosMultiphysics.FluidDynamicsApplication import FSStrategy as solver_strategy
    from KratosMultiphysics.FluidDynamicsApplication import StrategyLabel as strategy_label
    from KratosMultiphysics.RANSApplication import CoupledStrategyItem as coupled_strategy_item
else:
    raise Exception("Distributed run requires TrilinosApplication")

def CreateFractionalStepSolverSettings(
                    is_periodic,
                    model_part,
                    domain_size,
                    time_order,
                    use_slip_conditions,
                    move_mesh_flag,
                    reform_dofs_at_each_step,
                    communicator):
    args = []
    if (IsDistributedRun()):
        args.append(communicator)
    args.extend([model_part, domain_size, time_order, use_slip_conditions, move_mesh_flag, reform_dofs_at_each_step])

    if (is_periodic):
        args.append(KratosCFD.PATCH_INDEX)
        return periodic_solver_settings(*args)
    else:
        return solver_settings(*args)

class FractionalStepFormulation(Formulation):
    def __init__(self, model_part, settings):
        super(FractionalStepFormulation, self).__init__(model_part, settings)

        ##settings string in json format
        default_settings = Kratos.Parameters("""
        {
            "predictor_corrector": false,
            "maximum_velocity_iterations": 3,
            "maximum_pressure_iterations": 3,
            "velocity_tolerance": 1e-3,
            "pressure_tolerance": 1e-2,
            "dynamic_tau": 0.01,
            "oss_switch": 0,
            "echo_level": 0,
            "time_order": 2,
            "compute_reactions": false,
            "reform_dofs_at_each_step": false,
            "use_slip_conditions": true,
            "pressure_linear_solver_settings":  {
                "solver_type"                    : "amgcl",
                "max_iteration"                  : 200,
                "tolerance"                      : 1e-6,
                "provide_coordinates"            : false,
                "smoother_type"                  : "ilu0",
                "krylov_type"                    : "cg",
                "gmres_krylov_space_dimension"   : 100,
                "use_block_matrices_if_possible" : false,
                "coarsening_type"                : "aggregation",
                "scaling"                        : true,
                "verbosity"                      : 0
            },
            "velocity_linear_solver_settings": {
                "solver_type"                    : "amgcl",
                "max_iteration"                  : 200,
                "tolerance"                      : 1e-6,
                "provide_coordinates"            : false,
                "smoother_type"                  : "ilu0",
                "krylov_type"                    : "lgmres",
                "gmres_krylov_space_dimension"   : 100,
                "use_block_matrices_if_possible" : false,
                "coarsening_type"                : "aggregation",
                "scaling"                        : true,
                "verbosity"                      : 0
            },
            "update_processes_list": []
        }""")

        settings.ValidateAndAssignDefaults(default_settings)

        self.min_buffer_size = 3
        self.element_has_nodal_properties = True
        self.fractional_step_model_part = None

        ## Construct the linear solvers
        self.pressure_linear_solver = CreateLinearSolver(self.settings["pressure_linear_solver_settings"])
        self.velocity_linear_solver = CreateLinearSolver(self.settings["velocity_linear_solver_settings"])

        self.compute_reactions = self.settings["compute_reactions"].GetBool()

        Kratos.Logger.PrintInfo("FractionalStepFormulation", "Construction of FractionalStepFormulation finished.")

    def AddVariables(self):
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.DENSITY)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.PRESSURE)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.VELOCITY)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.ACCELERATION)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.MESH_VELOCITY)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.BODY_FORCE)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.NODAL_H)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.REACTION)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.REACTION_WATER_PRESSURE)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.NORMAL)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.Y_WALL)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.EXTERNAL_PRESSURE)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.VISCOSITY)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.FRACT_VEL)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.PRESSURE_OLD_IT)
        # The following are used for the calculation of projections
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.NODAL_AREA)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.PRESS_PROJ)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.CONV_PROJ)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.DIVPROJ)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(KratosCFD.Q_VALUE)

        if self.IsPeriodic():
            self.GetBaseModelPart().AddNodalSolutionStepVariable(KratosCFD.PATCH_INDEX)

        Kratos.Logger.PrintInfo("FractionalStepFormulation", "Added solution step variables.")

    def AddDofs(self):
        Kratos.VariableUtils().AddDof(Kratos.VELOCITY_X, Kratos.REACTION_X,self.GetBaseModelPart())
        Kratos.VariableUtils().AddDof(Kratos.VELOCITY_Y, Kratos.REACTION_Y,self.GetBaseModelPart())
        Kratos.VariableUtils().AddDof(Kratos.VELOCITY_Z, Kratos.REACTION_Z,self.GetBaseModelPart())
        Kratos.VariableUtils().AddDof(Kratos.PRESSURE, Kratos.REACTION_WATER_PRESSURE,self.GetBaseModelPart())

    def PrepareModelPart(self):
        self.domain_size = self.GetBaseModelPart().ProcessInfo[Kratos.DOMAIN_SIZE]
        condition_suffix = str(self.domain_size) + "D" + str(self.domain_size) + "N"
        element_suffix = str(self.domain_size) + "D" + str(self.domain_size + 1) + "N"

        # preparing velocity_pressure model part
        element_name = "FractionalStep" + element_suffix
        condition_name = "FSHighReKWall" + condition_suffix
        modelpart_name = self.GetBaseModelPart().Name + "_" + element_name + "_" + condition_name
        self.fractional_step_model_part =  CreateDuplicateModelPart(
                                                        self.GetBaseModelPart(),
                                                        modelpart_name,
                                                        element_name,
                                                        condition_name,
                                                        "")

    def Initialize(self):
        self.solver_settings = CreateFractionalStepSolverSettings(
                                            self.IsPeriodic(),
                                            self.fractional_step_model_part,
                                            self.domain_size,
                                            self.settings["time_order"].GetInt(),
                                            self.settings["use_slip_conditions"].GetBool(),
                                            self.GetMoveMeshFlag(),
                                            self.settings["reform_dofs_at_each_step"].GetBool(),
                                            self.GetCommunicator())

        self.solver_settings.SetEchoLevel(self.settings["echo_level"].GetInt())

        self.solver_settings.SetStrategy(strategy_label.Velocity,
                                         self.velocity_linear_solver,
                                         self.settings["velocity_tolerance"].GetDouble(),
                                         self.settings["maximum_velocity_iterations"].GetInt())

        self.solver_settings.SetStrategy(strategy_label.Pressure,
                                         self.pressure_linear_solver,
                                         self.settings["pressure_tolerance"].GetDouble(),
                                         self.settings["maximum_pressure_iterations"].GetInt())


        if self.IsPeriodic():
            self.solver = solver_strategy(self.fractional_step_model_part,
                                          self.solver_settings,
                                          self.settings["predictor_corrector"].GetBool(),
                                          KratosCFD.PATCH_INDEX)
        else:
            self.solver = solver_strategy(self.fractional_step_model_part,
                                          self.solver_settings,
                                          self.settings["predictor_corrector"].GetBool())

        self.GetBaseModelPart().ProcessInfo.SetValue(Kratos.DYNAMIC_TAU, self.settings["dynamic_tau"].GetDouble())
        self.GetBaseModelPart().ProcessInfo.SetValue(Kratos.OSS_SWITCH, self.settings["oss_switch"].GetInt())
        self.GetBaseModelPart().ProcessInfo.SetValue(Kratos.PRESSURE_COEFFICIENT, 0.5)

        Kratos.Logger.PrintInfo("FractionalStepFormulation", "Solver initialization finished.")

    def GetMinimumBufferSize(self):
        return self.min_buffer_size

    def Finalize(self):
        self.solver.Clear()

    def GetCoupledStrategyItems(self):
        current_strategy_item = coupled_strategy_item(
                                        self.solver,
                                        "FractionalStepVelocityPressure",
                                        self.settings["echo_level"].GetInt())
        return [["velocity_pressure_solver_settings", current_strategy_item]]

    def IsSolvingForSteadyState(self):
        return False