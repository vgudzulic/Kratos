# This class can be taken as an abstract template for derivation. It can be used
# for default passive behaviour.

from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.SwimmingDEMApplication as SDEM

class DerivativesRecoverer:
    def __init__(self, project_parameters, model_part):
        self.model_part = model_part
        self.cplusplus_recovery_tool = SDEM.DerivativeRecoveryTool3D(model_part, project_parameters)

class EmptyGradientRecoverer(DerivativesRecoverer):
    def __init__(self, project_parameters, model_part):
        DerivativesRecoverer.__init__(self, project_parameters, model_part)
    def RecoverGradientOfScalar(self, scalar_variable, gradient_variable):
        pass
    def RecoverGradientOfVector(self, vector_variable, gradient_variable_x, gradient_variable_y, gradient_variable_z):
        pass
    def RecoverGradientOfVelocity(self):
        pass
    def RecoverFluidFractionGradient(self):
        pass
    def RecoverPressureGradient(self):
        pass

class EmptyMaterialAccelerationRecoverer(DerivativesRecoverer):
    def __init__(self, project_parameters, model_part):
        pass
    def RecoverMaterialAcceleration(self):
        pass
    def RecoverMaterialAccelerationFromGradient(self):
        pass

class EmptyVorticityRecoverer(DerivativesRecoverer):
    def __init__(self, project_parameters, model_part):
        pass
    def RecoverVorticityFromGradient(self):
        pass
    def CalculateVorticityContributionOfTheGradientOfAComponent(self):
        pass

class EmptyLaplacianRecoverer(DerivativesRecoverer):
    def __init__(self, project_parameters, model_part):
        pass
    def RecoverVectorLaplacian(self, vector_variable, laplacian_variable):
        pass
    def RecoverVelocityLaplacian(self):
        pass

class GradientRecoverer(EmptyGradientRecoverer):
    def __init__(self, project_parameters, model_part):
        DerivativesRecoverer.__init__(self, project_parameters, model_part)
    def RecoverGradientOfVelocity(self):
        self.RecoverGradientOfVector(KM.VELOCITY, KM.VELOCITY_X_GRADIENT, KM.VELOCITY_Y_GRADIENT, KM.VELOCITY_Z_GRADIENT)
    def RecoverPressureGradient(self):
        self.RecoverGradientOfScalar(KM.PRESSURE, KM.PRESSURE_GRADIENT)
    def RecoverFluidFractionGradient(self):
        self.RecoverGradientOfScalar(KM.FLUID_FRACTION, KM.FLUID_FRACTION_GRADIENT)

class MaterialAccelerationRecoverer(GradientRecoverer, EmptyMaterialAccelerationRecoverer):
    def __init__(self, project_parameters, model_part):
        GradientRecoverer.__init__(self, project_parameters, model_part)
    def RecoverMaterialAcceleration(self):
        self.RecoverMaterialAccelerationFromGradient()
    def RecoverMaterialAccelerationFromGradient(self):
        self.cplusplus_recovery_tool.CalculateVectorMaterialDerivativeFromGradient(self.model_part,
                                                                                   KM.VELOCITY_X_GRADIENT,
                                                                                   KM.VELOCITY_Y_GRADIENT,
                                                                                   KM.VELOCITY_Z_GRADIENT,
                                                                                   KM.ACCELERATION,
                                                                                   KM.MATERIAL_ACCELERATION)

class VorticityRecoverer(GradientRecoverer, EmptyVorticityRecoverer):
    def __init__(self, project_parameters, model_part):
        GradientRecoverer.__init__(self, project_parameters, model_part)
    def RecoverVorticityFromGradient(self):
        self.cplusplus_recovery_tool.CalculateVorticityFromGradient(self.model_part,
                                                                    KM.VELOCITY_X_GRADIENT,
                                                                    KM.VELOCITY_Y_GRADIENT,
                                                                    KM.VELOCITY_Z_GRADIENT,
                                                                    KM.VORTICITY)

    def CalculateVorticityContributionOfTheGradientOfAComponent(self):
        self.cplusplus_recovery_tool.CalculateVorticityContributionOfTheGradientOfAComponent(self.model_part,
                                                                                             KM.VELOCITY_COMPONENT_GRADIENT,
                                                                                             KM.VORTICITY)

class LaplacianRecoverer(GradientRecoverer, EmptyLaplacianRecoverer):
    def __init__(self, project_parameters, model_part):
        GradientRecoverer.__init__(self, project_parameters, model_part)
    def RecoverVelocityLaplacian(self):
        self.RecoverVectorLaplacian(KM.VELOCITY, KM.VELOCITY_LAPLACIAN)
