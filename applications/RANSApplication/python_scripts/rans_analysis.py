from __future__ import absolute_import, division  #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from sys import argv

import KratosMultiphysics as Kratos
import KratosMultiphysics.RANSApplication as KratosRANS
from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis
from KratosMultiphysics.RANSApplication import RansVariableUtilities


class RANSAnalysis(FluidDynamicsAnalysis):
    '''Main script for fluid dynamics simulations using the navier_stokes family of python solvers.'''
    def __init__(self, model, parameters):
        super(RANSAnalysis, self).__init__(model, parameters)

    def KeepAdvancingSolutionLoop(self):
        velocity_convergence = False
        pressure_convergence = False
        if self._GetSolver()._TimeBufferIsInitialized():
            velocity_convergence = RansVariableUtilities.CalculateTransientVariableConvergence(
                self._GetSolver().GetComputingModelPart(), Kratos.VELOCITY,
                1e-3, 1e-5, 1)
            pressure_convergence = RansVariableUtilities.CalculateTransientVariableConvergence(
                self._GetSolver().GetComputingModelPart(), Kratos.PRESSURE,
                1e-3, 1e-5, 1)
            tke_convergence = RansVariableUtilities.CalculateTransientVariableConvergence(
                self._GetSolver().GetComputingModelPart(), KratosRANS.TURBULENT_KINETIC_ENERGY,
                1e-3, 1e-5, 1)
            epsilon_convergence = RansVariableUtilities.CalculateTransientVariableConvergence(
                self._GetSolver().GetComputingModelPart(), KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE,
                1e-3, 1e-5, 1)

            if (velocity_convergence and pressure_convergence and tke_convergence and epsilon_convergence):
                Kratos.Logger.PrintInfo(self.__class__.__name__,
                                        "Steady state solution reached.")
        return ((self.time < self.end_time)
                and not (velocity_convergence and pressure_convergence and tke_convergence and epsilon_convergence))

    def _GetSimulationName(self):
        return "RANS Analysis"


if __name__ == '__main__':
    if len(argv) > 2:
        err_msg = 'Too many input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '- With default parameter file (assumed to be called "ProjectParameters.json"):\n'
        err_msg += '    "python fluid_dynamics_analysis.py"\n'
        err_msg += '- With custom parameter file:\n'
        err_msg += '    "python fluid_dynamics_analysis.py <my-parameter-file>.json"\n'
        raise Exception(err_msg)

    if len(argv) == 2:  # ProjectParameters is being passed from outside
        parameter_file_name = argv[1]
    else:  # using default name
        parameter_file_name = "ProjectParameters.json"

    with open(parameter_file_name, 'r') as parameter_file:
        parameters = Kratos.Parameters(parameter_file.read())

    model = Kratos.Model()
    simulation = RANSAnalysis(model, parameters)
    simulation.Run()
