import KratosMultiphysics
import KratosMultiphysics.FluidTransportApplication as KratosFluidTransport

from math import pi
from sympy import Symbol, diff, lambdify, sin, cos, atan, exp
from sympy.vector import CoordSys3D, gradient, divergence #, laplacian
import numpy as np

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyScalarManufacturedConstraintFunctionProcess(Model, settings["Parameters"])

## All the python processes should be derived from "python_process"

class ApplyScalarManufacturedConstraintFunctionProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        self.model_part = Model[settings["model_part_name"].GetString()]

        # self.a1 = 2
        # self.a2 = 3
        k = 0.000001
        s = 20

        time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        for node in self.model_part.Nodes:
            f = 4.0*node.X*node.Y*k**(-0.5)*(1.0 - 2*node.X)*(1 - node.Y)*(16 - 16*node.X)*sin(pi*time)/(pi*(4*k**(-1.0)*(-(node.X - 0.5)**2 - (node.Y - 0.5)**2 + 0.0625)**2 + 1)) + 6.0*node.X*node.Y*k**(-0.5)*(1.0 - 2*node.Y)*(1 - node.Y)*(16 - 16*node.X)*sin(pi*time)/(pi*(4*k**(-1.0)*(-(node.X - 0.5)**2 - (node.Y - 0.5)**2 + 0.0625)**2 + 1)) + 16*node.X*node.Y*s*(1 - node.X)*(1 - node.Y)*(atan(2*k**(-0.5)*(-(node.X - 0.5)**2 - (node.Y - 0.5)**2 + 0.0625))/pi + 0.5)*sin(pi*time) + node.X*node.Y*pi*(1 - node.Y)*(16 - 16*node.X)*(atan(2*k**(-0.5)*(-(node.X - 0.5)**2 - (node.Y - 0.5)**2 + 0.0625))/pi + 0.5)*cos(pi*time) + 3.0*node.X*node.Y*(16*node.X - 16)*(atan(2*k**(-0.5)*(-(node.X - 0.5)**2 - (node.Y - 0.5)**2 + 0.0625))/pi + 0.5)*sin(pi*time) + 2.0*node.X*node.Y*(16*node.Y - 16)*(atan(2*k**(-0.5)*(-(node.X - 0.5)**2 - (node.Y - 0.5)**2 + 0.0625))/pi + 0.5)*sin(pi*time) + 3.0*node.X*(1 - node.Y)*(16 - 16*node.X)*(atan(2*k**(-0.5)*(-(node.X - 0.5)**2 - (node.Y - 0.5)**2 + 0.0625))/pi + 0.5)*sin(pi*time) + 2.0*node.Y*(1 - node.Y)*(16 - 16*node.X)*(atan(2*k**(-0.5)*(-(node.X - 0.5)**2 - (node.Y - 0.5)**2 + 0.0625))/pi + 0.5)*sin(pi*time) - k*(-8*node.X*node.Y*k**(-1.5)*(1.0 - 2*node.X)*(1 - node.Y)*(2.0 - 4*node.X)*(16 - 16*node.X)*(-(node.X - 0.5)**2 - (node.Y - 0.5)**2 + 0.0625)*sin(pi*time)/(pi*(4*k**(-1.0)*(-(node.X - 0.5)**2 - (node.Y - 0.5)**2 + 0.0625)**2 + 1)**2) - 8*node.X*node.Y*k**(-1.5)*(1.0 - 2*node.Y)*(1 - node.Y)*(2.0 - 4*node.Y)*(16 - 16*node.X)*(-(node.X - 0.5)**2 - (node.Y - 0.5)**2 + 0.0625)*sin(pi*time)/(pi*(4*k**(-1.0)*(-(node.X - 0.5)**2 - (node.Y - 0.5)**2 + 0.0625)**2 + 1)**2) - 32*node.X*node.Y*k**(-0.5)*(1.0 - 2*node.X)*(1 - node.Y)*sin(pi*time)/(pi*(4*k**(-1.0)*(-(node.X - 0.5)**2 - (node.Y - 0.5)**2 + 0.0625)**2 + 1)) + 2*node.X*node.Y*k**(-0.5)*(1.0 - 2*node.X)*(16*node.Y - 16)*sin(pi*time)/(pi*(4*k**(-1.0)*(-(node.X - 0.5)**2 - (node.Y - 0.5)**2 + 0.0625)**2 + 1)) - 2*node.X*node.Y*k**(-0.5)*(1.0 - 2*node.Y)*(16 - 16*node.X)*sin(pi*time)/(pi*(4*k**(-1.0)*(-(node.X - 0.5)**2 - (node.Y - 0.5)**2 + 0.0625)**2 + 1)) + 2*node.X*node.Y*k**(-0.5)*(1.0 - 2*node.Y)*(16*node.X - 16)*sin(pi*time)/(pi*(4*k**(-1.0)*(-(node.X - 0.5)**2 - (node.Y - 0.5)**2 + 0.0625)**2 + 1)) - 8*node.X*node.Y*k**(-0.5)*(1 - node.Y)*(16 - 16*node.X)*sin(pi*time)/(pi*(4*k**(-1.0)*(-(node.X - 0.5)**2 - (node.Y - 0.5)**2 + 0.0625)**2 + 1)) + 4*node.X*k**(-0.5)*(1.0 - 2*node.Y)*(1 - node.Y)*(16 - 16*node.X)*sin(pi*time)/(pi*(4*k**(-1.0)*(-(node.X - 0.5)**2 - (node.Y - 0.5)**2 + 0.0625)**2 + 1)) + 2*node.X*(16*node.X - 16)*(atan(2*k**(-0.5)*(-(node.X - 0.5)**2 - (node.Y - 0.5)**2 + 0.0625))/pi + 0.5)*sin(pi*time) + 4*node.Y*k**(-0.5)*(1.0 - 2*node.X)*(1 - node.Y)*(16 - 16*node.X)*sin(pi*time)/(pi*(4*k**(-1.0)*(-(node.X - 0.5)**2 - (node.Y - 0.5)**2 + 0.0625)**2 + 1)) + 2*node.Y*(16*node.Y - 16)*(atan(2*k**(-0.5)*(-(node.X - 0.5)**2 - (node.Y - 0.5)**2 + 0.0625))/pi + 0.5)*sin(pi*time))
            node.SetSolutionStepValue(KratosMultiphysics.HEAT_FLUX,f)
        # self.components_process_list = []

        # for node in model_part.Nodes:
        #     temperature = 5.0 / 8.0 * node.X + 3.0
        #     node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE,temperature)


    # def ExecuteInitialize(self):

    #     for component in self.components_process_list:
    #         component.ExecuteInitialize()

    # def ExecuteInitializeSolutionStep(self):

    #     for component in self.components_process_list:
    #         component.ExecuteInitializeSolutionStep()

    # this is a draft
# a Kratos Process shall be developet to apply the source term at the InitializeSolutionStep




    # def ExecuteInitializeSolutionStep(self):
    #     time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
    #     for node in self.model_part.Nodes:
    #         f = self.source_term(node.X, node.Y, node.Z, time)
    #         a = self.velocity(node.X, node.Y, node.Z, time)
    #         node.SetSolutionStepValue(KratosMultiphysics.HEAT_FLUX,f)
    #         # node.FastGetSolutionStepValue(HEAT_FLUX) = f
    #        # node.FastGetSolutionStepValue(VELOCITY) = a

    # # def ExecuteFinalizeSolutionStep(self):
    # #     time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
    # #     for node in self.model_part.Nodes:
    # #         u = self.analytic_phi(node.X, node.Y, node.Z, time)
    # #         node.FastGetSolutionStepValue(ANALYTIC_PHI) = u

    # def _GenerateAnalyticFields(self, settings, c, t):
    #     # a1 = settings['velocity'][0].GetDouble()
    #     # a2 = settings['velocity'][1].GetDouble()
    #     # k = settings['diffusivity'].GetDouble()
    #     # s = settings['absorption'].GetDouble()
    #     k = 0.000001
    #     s = 1000
    #     # a = a1*c.i + a2*c.j
    #     a =[2.0, 3.0]
    #     u = 16 * sin(pi * t) * c.x*(1 - c.x) * c.y*(1 - c.y) * (1/2 + (atan(2 * k ** (-1/2) * (0.25**2 - (c.x - 0.5)**2 - (c.y - 0.5)**2)))/pi)

    #     return a, u
