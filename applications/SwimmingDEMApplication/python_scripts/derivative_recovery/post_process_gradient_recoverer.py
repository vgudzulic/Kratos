from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.DEMApplication as DEM
import KratosMultiphysics.SwimmingDEMApplication as SDEM
from . import recoverer

class PostProcessGradientMaterialAccelerationRecoverer(recoverer.MaterialAccelerationRecoverer):
    def __init__(self, project_parameters, model_part):
        self.store_full_gradient = project_parameters["store_full_gradient_option"].GetBool()
        recoverer.MaterialAccelerationRecoverer.__init__(self, project_parameters, model_part)

    def RecoverMaterialAcceleration(self):
        self.cplusplus_recovery_tool.CalculateVectorMaterialDerivative(self.model_part,
                                                                       KM.VELOCITY,
                                                                       KM.ACCELERATION,
                                                                       KM.MATERIAL_ACCELERATION)
