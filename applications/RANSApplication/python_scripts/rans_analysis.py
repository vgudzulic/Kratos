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

    def RunSolutionLoop(self):
        """This function executes the solution loop of the AnalysisStage
        It can be overridden by derived classes
        """
        is_converged = False
        while (self.KeepAdvancingSolutionLoop() and not is_converged):
            self.time = self._GetSolver().AdvanceInTime(self.time)
            self.InitializeSolutionStep()
            self._GetSolver().Predict()
            is_converged = self._GetSolver().SolveSolutionStep()
            self.FinalizeSolutionStep()
            self.OutputSolutionStep()

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
