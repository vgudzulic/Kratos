import KratosMultiphysics
import KratosMultiphysics.RANSConstitutiveLawsApplication as KratosRANS

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyKEpsilonUtauProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ApplyKEpsilonUtauProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters( """
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "y_plus_calculation"      : {
                    "model_type"     : "logarithmic"
                },
                "constants"            : {},
                "is_initialization_process" : false
            }  """ )

        settings.ValidateAndAssignDefaults(default_parameters)
        self.settings = settings

        self.model_part = Model[settings["model_part_name"].GetString()]
        from rans_y_plus_model_factory import Factory

        self.y_plus_model = Factory(self.model_part, self.settings["y_plus_calculation"])
        self.settings.RemoveValue("y_plus_calculation")


        self.k_epsilon_utau_process = KratosRANS.RansKEpsilonEvaluationUtauProcess(Model, self.settings, self.y_plus_model)

    def Check(self):
        self.y_plus_model.Check()
        self.k_epsilon_utau_process.Check()

    def ExecuteInitializeSolutionStep(self):
        if (self.model_part.ProcessInfo[KratosRANS.IS_CO_SOLVING_PROCESS_ACTIVE]):
            self.k_epsilon_utau_process.Execute()
