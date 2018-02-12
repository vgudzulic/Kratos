from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("ALEApplication", "TrilinosApplication")

# Import applications
import KratosMultiphysics.TrilinosApplication as TrilinosApplication

# Other imports
from KratosMultiphysics.mpi import *
import os
import trilinos_mesh_solver_base


def CreateSolver(mesh_model_part, custom_settings):
    return TrilinosMeshSolverStructuralSimilarity(mesh_model_part, custom_settings)


class TrilinosMeshSolverStructuralSimilarity(trilinos_mesh_solver_base.TrilinosMeshSolverBase):
    def __init__(self, mesh_model_part, custom_settings):
        super(TrilinosMeshSolverStructuralSimilarity, self).__init__(mesh_model_part, custom_settings)
        mpi.world.barrier()
        if mpi.rank == 0:
            print("::[TrilinosMeshSolverStructuralSimilarity]:: Construction finished")

    #### Private functions ####

    def _create_mesh_motion_solver(self):
        linear_solver = self.get_linear_solver()
        communicator = self.get_communicator()
        time_order = self.settings["time_order"].GetInt()
        reform_dofs_each_step = self.settings["reform_dofs_each_step"].GetBool()
        compute_reactions = self.settings["compute_reactions"].GetBool()
        echo_level = self.settings["echo_level"].GetInt()
        solver = TrilinosApplication.TrilinosStructuralMeshMovingStrategy(
            communicator,
            self.mesh_model_part,
            linear_solver,
            time_order,
            reform_dofs_each_step,
            compute_reactions,
            echo_level)
        return solver