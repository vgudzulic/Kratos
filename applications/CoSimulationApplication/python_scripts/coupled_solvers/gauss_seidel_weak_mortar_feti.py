# Importing the base class
from KratosMultiphysics.CoSimulationApplication.coupled_solvers.gauss_seidel_weak import GaussSeidelWeakCoupledSolver
# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_coupled_solver import CoSimulationCoupledSolver


def Create(settings, models, solver_name):
    return GaussSeidelWeakMortarFETICoupledSolver(settings, models, solver_name)

class GaussSeidelWeakMortarFETICoupledSolver(GaussSeidelWeakCoupledSolver):
    def __init__(self, settings, models, solver_name):
        print('\n\n ============ GaussSeidelWeakMortarFETICoupledSolver ============\n\n')
        super(GaussSeidelWeakMortarFETICoupledSolver, self).__init__(settings, models, solver_name)
        self.is_got_system_matrices = False
        print(self.data_transfer_operators_dict)
        self.system_tangents = []

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()
        if self.is_got_system_matrices == False:
            for solver_name, solver in self.solver_wrappers.items():
                print('solver = ',)
                self.system_tangents.append(solver.get_mechanical_solution_strategy().GetSystemMatrix())
        print('\n\n ------- PRINTING SYSTEM MATRICES \n',self.system_tangents)


    def SolveSolutionStep(self):

        for coupling_op in self.coupling_operations_dict.values():
            coupling_op.InitializeCouplingIteration()

        for solver_name, solver in self.solver_wrappers.items():
            self._SynchronizeInputData(solver_name)
            solver.SolveSolutionStep()
            self._SynchronizeOutputData(solver_name)

        for coupling_op in self.coupling_operations_dict.values():
            coupling_op.FinalizeCouplingIteration()

        return True
