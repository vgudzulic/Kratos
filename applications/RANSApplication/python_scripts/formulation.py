from abc import ABC, abstractmethod


class Formulation(ABC):
    def __init__(self, base_computing_model_part, settings):
        self.settings = settings
        self.base_computing_model_part = base_computing_model_part
        self.communicator = None
        self.is_periodic = False
        self.move_mesh = False

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
    def GetMinimumBufferSize(self):
        return 0

    @abstractmethod
    def Finalize(self):
        pass

    @abstractmethod
    def IsSolvingForSteadyState(self):
        return False

    @abstractmethod
    def GetCoupledStrategyItems(self):
        return None

    def IsPeriodic(self):
        return self.is_periodic

    def SetIsPeriodic(self, value):
        self.is_periodic = value

    def SetCommunicator(self, communicator):
        self.communicator = communicator

    def SetMoveMeshFlag(self, value):
        self.move_mesh = value

    def GetMoveMeshFlag(self):
        return self.move_mesh

    def GetCommunicator(self):
        return self.communicator

    def GetBaseModelPart(self):
        return self.base_computing_model_part



