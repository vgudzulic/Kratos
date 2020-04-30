import KratosMultiphysics
import KratosMultiphysics.FluidTransportApplication as KratosFluidTransport

from math import pi
from sympy import Symbol, diff, lambdify, sin, atan, exp
from sympy.vector import CoordSysCartesian, gradient, divergence #, laplacian

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
        # self.k = 0.000001
        # self.s = 1000

        c = CoordSysCartesian('c')
        t = Symbol('t')

        # self.a = a1*c.i + a2*c.j
        # self.u = 16 * sympy.sin(pi * self.t) * c.x(1 - c.x) * c.y(1 - c.y) * (1/2 + (sympy.atan(2 * self.k ** (-1/2) * (0.25**2 - (c.x - 0.5)**2 - (c.y - 0.5)**2)))/pi)
        #  self.u = (c.x**2.0/2/a1 + k*c.x/a1**2 + (1.0/2/a1 + k/a1**2) * (exp(-a1/k) - exp(-a1/k*(1-c.x))) / (1 - exp(-a1/k)))
        a, u = self._GenerateAnalyticFields(settings, c, t)

        f = diff(u, t) - k*divergence(gradient(u, c), c) + a.dot(gradient(u, c)) + s*u

        # Converting the symbolic fields to numeric functions
        self.velocity = lambdify([c.x, c.y, c.z, t], [a.dot(c.i), a.dot(c.j), a.dot(c.k)], 'numpy')
        self.analytic_phi = lambdify([c.x, c.y, c.z, t], u, 'numpy')
        self.source_term = lambdify([c.x, c.y, c.z, t], f, 'numpy')

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
            f = self.source_term(node.X, node.Y, node.Z, time)
            a = self.velocity(node.X, node.Y, node.Z, time)
            node.FastGetSolutionStepValue(SOURCE_TERM) = f
            node.FastGetSolutionStepValue(VELOCITY) = a

    def ExecuteFinalizeSolutionStep(self):
        time = self.model_part.ProcessInfo[TIME]
        for node in self.model_part.Nodes:
            u = self.analytic_phi(node.X, node.Y, node.Z, time)
            node.FastGetSolutionStepValue(ANALYTIC_PHI) = u

    def _GenerateAnalyticFields(self, settings, c, t):
        a1 = settings['velocity'][0].GetDouble()
        a2 = settings['velocity'][1].GetDouble()
        k = settings['diffusivity'].GetDouble()
        s = settings['absorption'].GetDouble()

        a = a1*c.i + a2*c.j
        u = 16 * sin(pi * t) * c.x(1 - c.x) * c.y(1 - c.y) * (1/2 + (atan(2 * k ** (-1/2) * (0.25**2 - (c.x - 0.5)**2 - (c.y - 0.5)**2)))/pi)

        return a, u
