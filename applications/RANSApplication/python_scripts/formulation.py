from abc import ABC, abstractmethod


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

    @abstractmethod
    def AddVariables(self):
        pass

    @abstractmethod
    def AddDofs(self):
        pass

    @abstractmethod
    def PrepareModelPart(self):
        pass

    @abstractmethod
    def Initialize(self):
        pass

    @abstractmethod
    def IsPeriodic(self):
        pass

    def SetCommunicator(self, communicator):
        self.communicator = communicator

    def GetBaseModelPart(self):
        return self.base_model_part

    def GetStrategyList(self):
        return self.strategy_list
