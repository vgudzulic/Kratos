# Importing the Kratos Library
import KratosMultiphysics as KM

# CoSimulation imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
from KratosMultiphysics.CoSimulationApplication.coupling_interface_data import CouplingInterfaceData

def AllocateHistoricalVariablesFromCouplingData(data_settings, model, solver_name):
    '''This function retrieves the historical variables that are needed for the ModelParts from the
    specified CouplingInterfaceDatas and allocates them on the ModelParts
    Note that it can only be called after the (Main-)ModelParts are created
    '''
    data_settings = data_settings.Clone() # to not mess with the later validation
    for data_name, data_config in data_settings.items():
        print(data_name)
    stop
    for data in data_list:
        hist_var_dict = data.GetHistoricalVariableDict()
        for full_model_part_name, variable in hist_var_dict.items():
            main_model_part_name = full_model_part_name.split(".")[0]
            if not model.HasModelPart(main_model_part_name):
                raise Exception('ModelPart "{}" does not exist in solver "{}"!'.format(main_model_part_name, solver_name))
            main_model_part = model[main_model_part_name]
            if not main_model_part.HasNodalSolutionStepVariable(variable):
                cs_tools.cs_print_info("CoSimTools", 'Allocating historical variable "{}" in ModelPart "{}" for solver "{}"'.format(variable.Name(), main_model_part_name, solver_name))
                main_model_part.AddNodalSolutionStepVariable(variable)

def CreateMainModelPartsFromCouplingData(data_settings, model, solver_name):
    '''This function creates the Main-ModelParts that are used in the specified CouplingInterfaceDatas
    '''

    for data in data_settings:
        print(data.PrettyPrintJsonString())
        data.ValidateAndAssignDefaults(CouplingInterfaceData.GetDefaultParameters())

        main_model_part_name = data["model_part_name"].GetString().split(".")[0]
        if not model.HasModelPart(main_model_part_name):
            model.CreateModelPart(main_model_part_name)
            cs_tools.cs_print_info("CoSimTools", 'Created ModelPart "{}" for solver "{}"'.format(main_model_part_name, solver_name))
    stophere

def RecursiveCreateModelParts(model_part, model_part_name):
    '''This function creates a hierarchy of SubModelParts on a given ModelPart
    '''
    model_part_name, *sub_model_part_names = model_part_name.split(".")
    if not model_part.HasSubModelPart(model_part_name):
        cs_tools.cs_print_info("CoSimTools", 'Created "{}" as SubModelPart of "{}"'.format(model_part_name, model_part.Name))
        model_part = model_part.CreateSubModelPart(model_part_name)
    if len(sub_model_part_names) > 0:
        RecursiveCreateModelParts(model_part, ".".join(sub_model_part_names))
