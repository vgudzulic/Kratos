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
                    "squeeze_amplitude" : 2.0,
                    "x1_origin" : 0.5,
                    "x2_origin" : 0.5,
                    "omega" : 0.0
            }
            """
            )

        settings.ValidateAndAssignDefaults(default_settings)

        self.L = settings["length"].GetDouble()
        self.rho = settings["density"].GetDouble()
        self.nu = settings["viscosity"].GetDouble() / self.rho
        self.delta_alpha = settings["delta_alpha"].GetDouble()
        self.omega = settings["omega"].GetDouble()