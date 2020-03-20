import KratosMultiphysics
import numpy as np

from KratosMultiphysics.SwimmingDEMApplication import field_utilities

## Import base class file
from KratosMultiphysics.SwimmingDEMApplication.custom_body_force.manufactured_solution import ManufacturedSolution

def CreateManufacturedSolution(custom_settings):
    return CasasFluidFractionSolution(custom_settings)

class CasasFluidFractionSolution(ManufacturedSolution):
    def __init__(self, settings):

        default_settings = KratosMultiphysics.Parameters("""
            {
                    "velocity"    : 1.0,
                    "length"      : 1.0,
                    "viscosity"   : 0.1,
                    "density"     : 1.0,
                    "frequency"   : 1.0,
                    "damping"     : 1.0,
                    "period"      : 0.1,
                    "delta_alpha" : 0.25,
                    "max_squeeze_fraction" : 0.5
            }
            """
            )

        settings.ValidateAndAssignDefaults(default_settings)

        self.L = settings["length"].GetDouble()
        self.rho = settings["density"].GetDouble()
        self.nu = settings["viscosity"].GetDouble() / self.rho
        self.delta_alpha = settings["delta_alpha"].GetDouble()
        period = settings["period"].GetDouble()
        self.omega = 2 * np.pi / period

    def alpha(self, time, x1, x2, x3):
        return 1 - 16*self.delta_alpha*(-self.L/(2*(self.max_squeeze_fraction*np.sin(self.omega*time) + 1)) + x2)*(self.L/(2*(self.max_squeeze_fraction*np.sin(self.omega*time) + 1)) + x2)*(-self.L*(self.max_squeeze_fraction*np.sin(self.omega*time) + 1)/2 + x1)*(self.L*(self.max_squeeze_fraction*np.sin(self.omega*time) + 1)/2 + x1)/self.L**4

    def alpha1(self, time, x1, x2, x3):
        return -16*self.delta_alpha*(-self.L/(2*(self.max_squeeze_fraction*np.sin(self.omega*time) + 1)) + x2)*(self.L/(2*(self.max_squeeze_fraction*np.sin(self.omega*time) + 1)) + x2)*(-self.L*(self.max_squeeze_fraction*np.sin(self.omega*time) + 1)/2 + x1)/self.L**4 - 16*self.delta_alpha*(-self.L/(2*(self.max_squeeze_fraction*np.sin(self.omega*time) + 1)) + x2)*(self.L/(2*(self.max_squeeze_fraction*np.sin(self.omega*time) + 1)) + x2)*(self.L*(self.max_squeeze_fraction*np.sin(self.omega*time) + 1)/2 + x1)/self.L**4

    def alpha2(self, time, x1, x2, x3):
        return -16*self.delta_alpha*(-self.L/(2*(self.max_squeeze_fraction*np.sin(self.omega*time) + 1)) + x2)*(-self.L*(self.max_squeeze_fraction*np.sin(self.omega*time) + 1)/2 + x1)*(self.L*(self.max_squeeze_fraction*np.sin(self.omega*time) + 1)/2 + x1)/self.L**4 - 16*self.delta_alpha*(self.L/(2*(self.max_squeeze_fraction*np.sin(self.omega*time) + 1)) + x2)*(-self.L*(self.max_squeeze_fraction*np.sin(self.omega*time) + 1)/2 + x1)*(self.L*(self.max_squeeze_fraction*np.sin(self.omega*time) + 1)/2 + x1)/self.L**4

    def alpha3(self, time, x1, x2, x3):
        return 0.0

    def dalphat(self, time, x1, x2, x3):
        return -8*self.delta_alpha*self.max_squeeze_fraction*self.omega*(-self.L/(2*(self.max_squeeze_fraction*np.sin(self.omega*time) + 1)) + x2)*(self.L/(2*(self.max_squeeze_fraction*np.sin(self.omega*time) + 1)) + x2)*(-self.L*(self.max_squeeze_fraction*np.sin(self.omega*time) + 1)/2 + x1)*np.cos(self.omega*time)/self.L**3 + 8*self.delta_alpha*self.max_squeeze_fraction*self.omega*(-self.L/(2*(self.max_squeeze_fraction*np.sin(self.omega*time) + 1)) + x2)*(self.L/(2*(self.max_squeeze_fraction*np.sin(self.omega*time) + 1)) + x2)*(self.L*(self.max_squeeze_fraction*np.sin(self.omega*time) + 1)/2 + x1)*np.cos(self.omega*time)/self.L**3 + 8*self.delta_alpha*self.max_squeeze_fraction*self.omega*(-self.L/(2*(self.max_squeeze_fraction*np.sin(self.omega*time) + 1)) + x2)*(-self.L*(self.max_squeeze_fraction*np.sin(self.omega*time) + 1)/2 + x1)*(self.L*(self.max_squeeze_fraction*np.sin(self.omega*time) + 1)/2 + x1)*np.cos(self.omega*time)/(self.L**3*(self.max_squeeze_fraction*np.sin(self.omega*time) + 1)**2) - 8*self.delta_alpha*self.max_squeeze_fraction*self.omega*(self.L/(2*(self.max_squeeze_fraction*np.sin(self.omega*time) + 1)) + x2)*(-self.L*(self.max_squeeze_fraction*np.sin(self.omega*time) + 1)/2 + x1)*(self.L*(self.max_squeeze_fraction*np.sin(self.omega*time) + 1)/2 + x1)*np.cos(self.omega*time)/(self.L**3*(self.max_squeeze_fraction*np.sin(self.omega*time) + 1)**2)