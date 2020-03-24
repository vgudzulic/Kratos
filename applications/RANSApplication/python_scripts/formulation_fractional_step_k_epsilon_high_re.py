# imports from kratos core
import KratosMultiphysics as Kratos
from KratosMultiphysics.process_factory import KratosProcessFactory
from KratosMultiphysics import python_linear_solver_factory as linear_solver_factory

# imports from fluid dynamics application
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import KratosMultiphysics.RANSApplication as KratosRANS

# imports from rans application
import KratosMultiphysics.RANSApplication as KratosRANS
from KratosMultiphysics.RANSApplication.formulation import Formulation
from KratosMultiphysics.RANSApplication.model_part_factory import CreateDuplicateModelPart
from KratosMultiphysics.RANSApplication.strategy_factory import CreateStrategy
from KratosMultiphysics.RANSApplication import RansCalculationUtilities


class FractionalStepKEpsilonHighRe(Formulation):
    def __init__(self, model, settings, scheme_settings):
        super(FractionalStepKEpsilonHighRe, self).__init__(model, settings, scheme_settings)
        defaults = Kratos.Parameters(r'''{
            "constants":
            {
                "wall_smoothness_beta"    : 5.2,
                "von_karman"              : 0.41,
                "c_mu"                    : 0.09,
                "c1"                      : 1.44,
                "c2"                      : 1.92,
                "sigma_k"                 : 1.0,
                "sigma_epsilon"           : 1.3,
                "dynamic_tau"             : 0.0,
                "oss_switch"              : 0
            },
            "velocity_pressure_solver_settings": {
                "predictor_corrector": false,
                "maximum_velocity_iterations": 3,
                "maximum_pressure_iterations": 3,
                "velocity_tolerance": 1e-3,
                "pressure_tolerance": 1e-2,
                "dynamic_tau": 0.01,
                "oss_switch": 0,
                "echo_level": 0,
                "consider_periodic_conditions": false,
                "time_order": 2,
                "compute_reactions": false,
                "reform_dofs_at_each_step": false,
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
                "move_mesh_flag": false,
                "use_slip_conditions": true
            },
            "turbulent_kinetic_energy_solver_settings": {
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
            "turbulent_kinetic_energy_solver_settings",
            "turbulent_energy_dissipation_rate_solver_settings"
            ]

        self.variable_list = [
            [KratosRANS.TURBULENT_KINETIC_ENERGY, KratosRANS.TURBULENT_KINETIC_ENERGY_RATE, KratosRANS.RANS_AUXILIARY_VARIABLE_1],
            [KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE, KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE_2, KratosRANS.RANS_AUXILIARY_VARIABLE_2]
        ]

    def AddVariables(self):
        self.base_model_part.AddNodalSolutionStepVariable(Kratos.DENSITY)
        self.base_model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        self.base_model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        self.base_model_part.AddNodalSolutionStepVariable(Kratos.ACCELERATION)
        self.base_model_part.AddNodalSolutionStepVariable(Kratos.MESH_VELOCITY)
        self.base_model_part.AddNodalSolutionStepVariable(Kratos.BODY_FORCE)
        self.base_model_part.AddNodalSolutionStepVariable(Kratos.NODAL_H)
        self.base_model_part.AddNodalSolutionStepVariable(Kratos.REACTION)
        self.base_model_part.AddNodalSolutionStepVariable(Kratos.REACTION_WATER_PRESSURE)
        self.base_model_part.AddNodalSolutionStepVariable(Kratos.NORMAL)
        self.base_model_part.AddNodalSolutionStepVariable(Kratos.Y_WALL)
        self.base_model_part.AddNodalSolutionStepVariable(Kratos.EXTERNAL_PRESSURE)
        self.base_model_part.AddNodalSolutionStepVariable(Kratos.VISCOSITY)
        self.base_model_part.AddNodalSolutionStepVariable(Kratos.FRACT_VEL)
        self.base_model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE_OLD_IT)

        # RANS specific variables
        self.base_model_part.AddNodalSolutionStepVariable(Kratos.RELAXED_ACCELERATION)
        self.base_model_part.AddNodalSolutionStepVariable(Kratos.KINEMATIC_VISCOSITY)
        self.base_model_part.AddNodalSolutionStepVariable(Kratos.TURBULENT_VISCOSITY)
        self.base_model_part.AddNodalSolutionStepVariable(KratosRANS.TURBULENT_KINETIC_ENERGY)
        self.base_model_part.AddNodalSolutionStepVariable(KratosRANS.TURBULENT_KINETIC_ENERGY_RATE)
        self.base_model_part.AddNodalSolutionStepVariable(KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE)
        self.base_model_part.AddNodalSolutionStepVariable(KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE_2)
        self.base_model_part.AddNodalSolutionStepVariable(KratosRANS.RANS_AUXILIARY_VARIABLE_1)
        self.base_model_part.AddNodalSolutionStepVariable(KratosRANS.RANS_AUXILIARY_VARIABLE_2)
        # The following are used for the calculation of projections
        self.base_model_part.AddNodalSolutionStepVariable(Kratos.NODAL_AREA)
        self.base_model_part.AddNodalSolutionStepVariable(Kratos.PRESS_PROJ)
        self.base_model_part.AddNodalSolutionStepVariable(Kratos.CONV_PROJ)
        self.base_model_part.AddNodalSolutionStepVariable(Kratos.DIVPROJ)

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

        # preparing velocity_pressure model part
        element_name = "FractionalStep" + element_suffix
        condition_name = "FSHighReKWall" + condition_suffix
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
        self.base_model_part.ProcessInfo[Kratos.DYNAMIC_TAU] = constants["dynamic_tau"].GetDouble()
        self.base_model_part.ProcessInfo[Kratos.OSS_SWITCH] = constants["dynamic_tau"].GetInt()
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
        if self.scheme_settings["scheme_type"].GetString() != "steady":
            raise Exception("Currently only steady schemes are supported")

        velocity_pressure_solver_settings = self.settings["velocity_pressure_solver_settings"]
        pressure_linear_solver = linear_solver_factory.ConstructSolver(velocity_pressure_solver_settings["pressure_linear_solver_settings"])
        velocity_linear_solver = linear_solver_factory.ConstructSolver(velocity_pressure_solver_settings["velocity_linear_solver_settings"])

        computing_model_part = self.model_part_list[0]
        #TODO: next part would be much cleaner if we passed directly the parameters to the c++
        fractionalstep_solver_settings = KratosCFD.FractionalStepSettings(computing_model_part,
                                                                            computing_model_part.ProcessInfo[Kratos.DOMAIN_SIZE],
                                                                            velocity_pressure_solver_settings["time_order"].GetInt(),
                                                                            velocity_pressure_solver_settings["use_slip_conditions"].GetBool(),
                                                                            velocity_pressure_solver_settings["move_mesh_flag"].GetBool(),
                                                                            velocity_pressure_solver_settings["reform_dofs_at_each_step"].GetBool())

        fractionalstep_solver_settings.SetEchoLevel(velocity_pressure_solver_settings["echo_level"].GetInt())

        fractionalstep_solver_settings.SetStrategy(KratosCFD.StrategyLabel.Velocity,
                                                   velocity_linear_solver,
                                                   velocity_pressure_solver_settings["velocity_tolerance"].GetDouble(),
                                                   velocity_pressure_solver_settings["maximum_velocity_iterations"].GetInt())

        fractionalstep_solver_settings.SetStrategy(KratosCFD.StrategyLabel.Pressure,
                                                   pressure_linear_solver,
                                                   velocity_pressure_solver_settings["pressure_tolerance"].GetDouble(),
                                                   velocity_pressure_solver_settings["maximum_pressure_iterations"].GetInt())


        solver_strategy = KratosCFD.FSStrategy(computing_model_part,
                                                fractionalstep_solver_settings,
                                                velocity_pressure_solver_settings["predictor_corrector"].GetBool())

        solver_strategy_dict = {}
        solver_strategy_dict["strategy"] = solver_strategy
        solver_strategy_dict["strategy_name"] = "VELOCITY_PRESSURE"
        solver_strategy_dict["update_processes_list"] = []
        self.strategy_list.append(solver_strategy_dict)

        default_strategy_settings = Kratos.Parameters(r'''{
            "update_processes": [],
            "strategy_settings": {}
        }''')

        for i, solver_name in enumerate(self.solver_names_list):
            settings = self.settings[solver_name]
            settings.ValidateAndAssignDefaults(default_strategy_settings)
            solver_strategy = {}
            solver_strategy["strategy"] = CreateStrategy(settings["strategy_settings"],
                                                         self.scheme_settings,
                                                         self.model_part_list[i+1],
                                                         self.base_model_part,
                                                         self.communicator,
                                                         self.variable_list[i][0],
                                                         self.variable_list[i][1],
                                                         self.variable_list[i][2])
            solver_strategy["strategy_name"] = self.variable_list[i][0].Name()

            factory = KratosProcessFactory(self.base_model_part.GetModel())
            solver_strategy["update_processes_list"] = factory.ConstructListOfProcesses(settings["update_processes"])

            self.strategy_list.append(solver_strategy)

    def IsPeriodic(self):
        is_periodic = False
        for solver_name in self.solver_names_list[1:]:
            if not is_periodic: is_periodic = self.settings[solver_name]["strategy_settings"]["is_periodic"].GetBool()
