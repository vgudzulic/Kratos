from __future__ import absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics.RANSApplication.coupled_rans_solver import CoupledRANSSolver as coupled_rans_solver


def CreateSolver(solver_name, model, settings):
    if (solver_name == "CoupledRANS"):
        return coupled_rans_solver(model, settings)
    else:
        raise Exception("Unsupported rans solver_name [ solver_name = \"" +
                        solver_name + "\" ]")
