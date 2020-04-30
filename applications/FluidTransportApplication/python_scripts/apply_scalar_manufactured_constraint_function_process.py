import KratosMultiphysics
import KratosMultiphysics.FluidTransportApplication as KratosFluidTransport
import math
import numpy as np
from sympy import Symbol, diff
from sympy.vector import CoordSysCartesian, gradient, divergence #, laplacian
from sympy.functions.elementary.exponential import exp

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyScalarManufacturedConstraintFunctionProcess(Model, settings["Parameters"])

## All the python processes should be derived from "python_process"

class ApplyScalarManufacturedConstraintFunctionProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        model_part = Model[settings["model_part_name"].GetString()]

        self.a1 = 2
        self.a2 = 3
        self.k = 0.000001
        self.s = 1000

        self.c = CoordSysCartesian('c')
        self.t = Symbol('t')

        self.a = a1*c.i + a2*c.j
        self.u = 16 * sympy.sin(math.pi * self.t) * c.x(1 - c.x) * c.y(1 - c.y) * (1/2 + (sympy.atan(2 * self.k ** (-1/2) * (0.25**2 - (c.x - 0.5)**2 - (c.y - 0.5)**2)))/math.pi)
        #  self.u = (c.x**2.0/2/a1 + k*c.x/a1**2 + (1.0/2/a1 + k/a1**2) * (exp(-a1/k) - exp(-a1/k*(1-c.x))) / (1 - exp(-a1/k)))

        self.f = diff(u, t) - k*divergence(gradient(u, c), c) + a.dot(gradient(u, c)) + s*u


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




    def ExecuteInitializeSolutionStep(self):
        time = self.model_part.ProcessInfo[TIME]
        for node in self.model_part.Nodes:
            f_analitica = self.f.eval('c.x'=node.X,...)
            node.FastGetSolutionStepValue(SOURCE_TERM) = f_analitica