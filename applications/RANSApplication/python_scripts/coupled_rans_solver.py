from __future__ import absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
from KratosMultiphysics.process_factory import KratosProcessFactory
from KratosMultiphysics.python_solver import PythonSolver

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import KratosMultiphysics.RANSApplication as KratosRANS

# Import application specific modules
from KratosMultiphysics.RANSApplication.formulation_factory import CreateFormulation
from KratosMultiphysics.RANSApplication.strategy_factory import CreateStrategy
from KratosMultiphysics.RANSApplication.incompressible_potential_flow_solver import IncompressiblePotentialFlowSolver
from KratosMultiphysics.FluidDynamicsApplication.check_and_prepare_model_process_fluid import CheckAndPrepareModelProcess


def CreateSolver(model, custom_settings):
    return CoupledRANSSolver(model, custom_settings)


class CoupledRANSSolver(PythonSolver):
    @classmethod
    def GetDefaultSettings(cls):

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "CoupledRANS",
            "model_part_name": "FluidModelPart",
            "domain_size": -1,
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name",
                "reorder": false
            },
            "material_import_settings": {
                "materials_filename": ""
            },
            "echo_level": 0,
            "volume_model_part_name": "volume_model_part",
            "skin_parts"   : [""],
            "no_skin_parts": [""],
            "assign_neighbour_elements_to_conditions": true,
            "time_stepping": {
                "automatic_time_step" : false,
                "CFL_number"          : 1,
                "minimum_delta_time"  : 1e-4,
                "maximum_delta_time"  : 0.01,
                "time_step"           : 0.0
            },
            "scheme_settings": {
                "scheme_type": "steady",
                "alpha_bossak": -0.3
            },
            "coupling_settings": {
                "formulation_type"    : "segregated_vms_k_epsilon_high_re",
                "max_iterations"      : 10,
                "formulation_settings": {}
            },
            "potential_flow_initialization": {}
        }""")

        default_settings.AddMissingParameters(super(CoupledRANSSolver, cls).GetDefaultSettings())
        return default_settings

    def __init__(self, model, custom_settings):
        self._validate_settings_in_baseclass=True # To be removed eventually
        super(CoupledRANSSolver, self).__init__(model,custom_settings)

        # Either retrieve the model part from the model or create a new one
        model_part_name = self.settings["model_part_name"].GetString()
        if model_part_name == "":
            raise Exception('Please provide the model part name as the "model_part_name" (string) parameter!')

        if self.model.HasModelPart(model_part_name):
            self.main_model_part = self.model.GetModelPart(model_part_name)
        else:
            self.main_model_part = self.model.CreateModelPart(model_part_name)

        domain_size = self.settings["domain_size"].GetInt()
        if domain_size == -1:
            raise Exception('Please provide the domain size as the "domain_size" (int) parameter!')
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, domain_size)

        formulation_name = self.settings["coupling_settings"]["formulation_type"].GetString()
        self.formulation = CreateFormulation(self.main_model_part,
                                             formulation_name,
                                             self.settings["coupling_settings"]["formulation_settings"],
                                             self.settings["scheme_settings"])

        if (not self.settings["potential_flow_initialization"].IsEquivalentTo(KratosMultiphysics.Parameters("{}"))):
            self.potential_flow_solver = IncompressiblePotentialFlowSolver(self.model, self.settings["potential_flow_initialization"])

        scheme_type = self.settings["scheme_settings"]["scheme_type"].GetString()
        if scheme_type == "bossak":
            self.min_buffer_size = 2
        elif scheme_type == "bdf3":
            self.min_buffer_size = 3
        elif scheme_type == "steady":
            self.min_buffer_size = 3
        else:
            msg  = "Unknown time_scheme option found in project parameters:\n"
            msg += "\"" + scheme_type + "\"\n"
            msg += "Accepted values are \"bossak\", or \"steady\".\n"
            raise Exception(msg)

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Solver construction finished.")

    def AddVariables(self):
        self.formulation.AddVariables()

        if (hasattr(self, "potential_flow_solver")):
            self.potential_flow_solver.fluid_model_part = self.main_model_part
            self.potential_flow_solver.AddVariables()

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Solver variables added correctly.")

    def AddDofs(self):
        self.formulation.AddDofs()

        if (hasattr(self, "potential_flow_solver")):
            self.potential_flow_solver.AddDofs()

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Solver dofs added correctly.")

    def ImportModelPart(self):
        # we can use the default implementation in the base class
        self._ImportModelPart(self.main_model_part, self.settings["model_import_settings"])

    def PrepareModelPart(self):
        if not self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
            ## Set fluid properties from materials json file
            materials_imported = self._SetPhysicalProperties()
            if not materials_imported:
                KratosMultiphysics.Logger.PrintWarning(self.__class__.__name__, "Material properties have not been imported. Check \'material_import_settings\' in your ProjectParameters.json.")
            ## Executes the check and prepare model process
            self._ExecuteCheckAndPrepare()
            ## Set buffer size
            self.main_model_part.SetBufferSize(self.min_buffer_size)

        if (hasattr(self, "potential_flow_solver")):
            self.potential_flow_solver.PrepareModelPart()

        self.formulation.PrepareModelPart()

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Model reading finished.")

    def ExportModelPart(self):
        ## Model part writing
        name_out_file = self.settings["model_import_settings"]["input_filename"].GetString()+".out"
        KratosMultiphysics.ModelPartIO(name_out_file, KratosMultiphysics.IO.WRITE).WriteModelPart(self.main_model_part)

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Model export finished.")

    def GetMinimumBufferSize(self):
        return self.min_buffer_size

    def Initialize(self):
        self.computing_model_part = self.GetComputingModelPart()

        # If needed, create the estimate time step utility
        if (self.settings["time_stepping"]["automatic_time_step"].GetBool()):
            self.EstimateDeltaTimeUtility = self._GetAutomaticTimeSteppingUtility()

        self.solver = KratosRANS.SegregatedStrategy(self.computing_model_part,
                                                    self.settings["coupling_settings"]["max_iterations"].GetInt())

        self.solver.SetEchoLevel(self.settings["echo_level"].GetInt())

        self.formulation.Initialize()

        # adding solving strategies
        for strategy_settings in self.formulation.GetStrategyList():
            strategy = strategy_settings["strategy"]
            strategy_name = strategy_settings["strategy_name"]

            self.solver.AddStrategy(strategy, strategy_name)
            for process in strategy_settings["update_processes_list"]:
                self.solver.AddAuxiliaryProcess(process, strategy_name)

        self.solver.Initialize()

        if (hasattr(self, "potential_flow_solver")):
            self.potential_flow_solver.Initialize()

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Solver initialization finished.")

    def AdvanceInTime(self, current_time):
        dt = self._ComputeDeltaTime()
        new_time = current_time + dt

        self.main_model_part.CloneTimeStep(new_time)
        self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] += 1

        return new_time

    def InitializeSolutionStep(self):
        if self._TimeBufferIsInitialized():
            self.solver.InitializeSolutionStep()

        if (hasattr(self, "potential_flow_solver")):
            self.potential_flow_solver.InitializeSolutionStep()

    def SolveSolutionStep(self):
        if self._TimeBufferIsInitialized():
            is_converged = self.solver.SolveSolutionStep()
            if not is_converged:
                msg  = "Fluid solver did not converge for step " + str(self.main_model_part.ProcessInfo[KratosMultiphysics.STEP]) + "\n"
                msg += "corresponding to time " + str(self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]) + "\n"
                KratosMultiphysics.Logger.PrintWarning(self.__class__.__name__, msg)
            return is_converged
        else:
            return True

    def FinalizeSolutionStep(self):
        if self._TimeBufferIsInitialized():
            self.solver.FinalizeSolutionStep()

    def Check(self):
        self.solver.Check()

        if (hasattr(self, "potential_flow_solver")):
            self.potential_flow_solver.Check()

    def Clear(self):
        self.solver.Clear()

        # if (hasattr(self, "potential_flow_solver")):
        #     self.potential_flow_solver.Clear()

    def GetComputingModelPart(self):
        if not self.main_model_part.HasSubModelPart("fluid_computational_model_part"):
            raise Exception("The ComputingModelPart was not created yet!")
        return self.main_model_part.GetSubModelPart("fluid_computational_model_part")

    def _TimeBufferIsInitialized(self):
        # We always have one extra old step (step 0, read from input)
        return self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] + 1 >= self.GetMinimumBufferSize()

    def _ComputeDeltaTime(self):
        # Automatic time step computation according to user defined CFL number
        if (self.settings["time_stepping"]["automatic_time_step"].GetBool()):
            delta_time = self.EstimateDeltaTimeUtility.EstimateDt()
        # User-defined delta time
        else:
            delta_time = self.settings["time_stepping"]["time_step"].GetDouble()

        return delta_time

    def _GetAutomaticTimeSteppingUtility(self):
        if (self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2):
            EstimateDeltaTimeUtility = KratosCFD.EstimateDtUtility2D(self.GetComputingModelPart(),
                                                                     self.settings["time_stepping"])
        else:
            EstimateDeltaTimeUtility = KratosCFD.EstimateDtUtility3D(self.GetComputingModelPart(),
                                                                     self.settings["time_stepping"])
        return EstimateDeltaTimeUtility

    def _ExecuteCheckAndPrepare(self):
        ## Check that the input read has the shape we like
        prepare_model_part_settings = KratosMultiphysics.Parameters("{}")
        prepare_model_part_settings.AddValue("volume_model_part_name", self.settings["volume_model_part_name"])
        prepare_model_part_settings.AddValue("skin_parts", self.settings["skin_parts"])
        if (self.settings.Has("assign_neighbour_elements_to_conditions")):
            prepare_model_part_settings.AddValue("assign_neighbour_elements_to_conditions",
                                                 self.settings["assign_neighbour_elements_to_conditions"])
        else:
            warn_msg = "\"assign_neighbour_elements_to_conditions\" should be added to defaults of " + self.__class__.__name__
            KratosMultiphysics.Logger.PrintWarning(
                '\n\x1b[1;31mDEPRECATION-WARNING\x1b[0m', warn_msg)

        CheckAndPrepareModelProcess(self.main_model_part, prepare_model_part_settings).Execute()

    def _SetPhysicalProperties(self):
        # Check if the fluid properties are provided using a .json file
        materials_filename = self.settings["material_import_settings"]["materials_filename"].GetString()
        if (materials_filename != ""):
            # Add constitutive laws and material properties from json file to model parts.
            material_settings = KratosMultiphysics.Parameters("""{"Parameters": {"materials_filename": ""}} """)
            material_settings["Parameters"]["materials_filename"].SetString(materials_filename)
            KratosMultiphysics.ReadMaterialsUtility(material_settings, self.model)
            materials_imported = True
        else:
            materials_imported = False

        return materials_imported