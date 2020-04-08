import KratosMultiphysics as Kratos
import KratosMultiphysics.RANSApplication as KratosRANS

from KratosMultiphysics import IsDistributedRun
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable

if (IsDistributedRun()
        and CheckIfApplicationsAvailable("TrilinosApplication")):
    from KratosMultiphysics.TrilinosApplication import trilinos_linear_solver_factory as linear_solver_factory
    from KratosMultiphysics.RANSApplication.TrilinosExtension import MPIGenericResidualBasedSimpleSteadyScalarScheme as steady_scheme
    from KratosMultiphysics.RANSApplication.TrilinosExtension import MPIResidualBasedSimpleSteadyVelocityScheme as steady_velocity_scheme
    from KratosMultiphysics.RANSApplication.TrilinosExtension import MPIGenericResidualBasedBossakVelocityDynamicScalarScheme as bossak_scheme
    from KratosMultiphysics.RANSApplication.TrilinosExtension import MPIGenericScalarConvergenceCriteria as scalar_convergence_criteria
    from KratosMultiphysics.TrilinosApplication import TrilinosResidualCriteria as residual_criteria
    from KratosMultiphysics.TrilinosApplication import TrilinosNewtonRaphsonStrategy as newton_raphson_strategy
    from KratosMultiphysics.RANSApplication.block_builder_and_solvers import TrilinosPeriodicBlockBuilderAndSolver as periodic_block_builder_and_solver
    from KratosMultiphysics.RANSApplication.block_builder_and_solvers import TrilinosBlockBuilderAndSolver as block_builder_and_solver
elif (not IsDistributedRun()):
    from KratosMultiphysics import python_linear_solver_factory as linear_solver_factory
    from KratosMultiphysics.RANSApplication import AlgebraicFluxCorrectedScalarSteadyScheme as steady_scheme
    from KratosMultiphysics.RANSApplication import ResidualBasedSimpleSteadyVelocityScheme as steady_velocity_scheme
    from KratosMultiphysics.RANSApplication import GenericResidualBasedBossakVelocityDynamicScalarScheme as bossak_scheme
    from KratosMultiphysics.RANSApplication import GenericScalarConvergenceCriteria as scalar_convergence_criteria
    from Kratos import ResidualCriteria as residual_criteria
    from Kratos import ResidualBasedNewtonRaphsonStrategy as newton_raphson_strategy
    from KratosMultiphysics.RANSApplication.block_builder_and_solvers import PeriodicBlockBuilderAndSolver as periodic_block_builder_and_solver
    from KratosMultiphysics.RANSApplication.block_builder_and_solvers import BlockBuilderAndSolver as block_builder_and_solver
else:
    raise Exception("Distributed run requires TrilinosApplication")


def CreateStrategy(solver_settings, scheme_settings, model_part,
                   base_model_part, communicator, variable,
                   variable_rate, relaxed_variable_rate):
    default_solver_settings = Kratos.Parameters(r'''{
            "is_periodic"           : false,
            "relative_tolerance"    : 1e-3,
            "absolute_tolerance"    : 1e-5,
            "max_iterations"        : 200,
            "relaxation_factor"     : 0.5,
            "echo_level"            : 0,
            "linear_solver_settings": {
                "solver_type"  : "amgcl"
            },
            "reform_dofs_at_each_step": true,
            "move_mesh_strategy": 0,
            "move_mesh_flag": false,
            "compute_reactions": false
    }''')

    default_scheme_settings = Kratos.Parameters(r'''{
        "scheme_type": "bossak",
        "scheme_settings": {}
    }''')

    solver_settings.ValidateAndAssignDefaults(default_solver_settings)
    scheme_settings.ValidateAndAssignDefaults(default_scheme_settings)

    linear_solver = linear_solver_factory.ConstructSolver(
        solver_settings["linear_solver_settings"])

    is_periodic = solver_settings["is_periodic"].GetBool()

    if is_periodic:
        InitializePeriodicConditions(model_part, base_model_part,
                                     variable)

    # TODO:
    if is_periodic and IsDistributedRun():
        msg = "\nCurrently periodic conditions in mpi is not supported due to following reasons:\n\n"
        msg += "    1. TrilinosResidualCriteria [ConvergenceCriterian]\n"
        msg += "PeriodicConditions duplicates one patch's equation ids to the counter patch. "
        msg += "The node and its corresponding dof might not fall in to the same partition raising an error in convergence calculation.\n\n"
        msg += "    2. ConnectivityPreserveModeller\n"
        msg += "Currently connectivity preserve modeller replaces all the conditions in an mdpa with given new condition. "
        msg += "This modeller is used to create modelparts having k-epsilon elements and conditions while sharing the same nodes as in VMS solution. "
        msg += "In the case of MPI, it is essential to have the PeriodicConditions in the mdpa file in order to properly distribute nodes to partitions using MetisApplication. "
        msg += "But if this is the case, PeriodicConditions also will be replaced by k-epsilon specific conditions casuing a segmentation fault.\n"
        msg += "    3. TrilinosBlockBuilderAndSolverPeriodic\n"
        msg += "In the case of MPI periodic in 2D, problem uses TrilinosBlockBuilderAndSolverPeriodic block builder and solver, which identifies "
        msg += "periodic conditions by number of nodes in the condition. So, In 2D all wall conditions and PeriodicConditions have only 2 nodes, all will be "
        msg += "considered as PeriodicConditions and will make the global assembly accordingly which is wrong."
        msg += "Therefore this error msg is printed in order to avoid confusion."
        raise Exception(msg)

    builder_and_solver = CreateBuilderAndSolver(linear_solver, is_periodic, communicator)

    if (scheme_settings["scheme_type"].GetString() == "transient"):
        convergence_criteria_type = scalar_convergence_criteria
    elif (scheme_settings["scheme_type"].GetString() == "steady"):
        convergence_criteria_type = residual_criteria

    convergence_criteria = convergence_criteria_type(
        solver_settings["relative_tolerance"].GetDouble(),
        solver_settings["absolute_tolerance"].GetDouble())

    if (scheme_settings["scheme_type"].GetString() == "transient"):
        time_scheme = bossak_scheme(scheme_settings["alpha_bossak"].GetDouble(),
                                     solver_settings["relaxation_factor"].GetDouble(),
                                     variable,
                                     variable_rate,
                                     relaxed_variable_rate)
    elif (scheme_settings["scheme_type"].GetString() == "steady"):
        base_model_part.ProcessInfo[Kratos.BOSSAK_ALPHA] = 0.0
        base_model_part.ProcessInfo[KratosRANS.IS_CO_SOLVING_PROCESS_ACTIVE] = True
        time_scheme = steady_scheme(solver_settings["relaxation_factor"].GetDouble())
    else:
        raise Exception("Unknown scheme_type = \"" +
                        scheme_settings["scheme_type"] + "\"")

    strategy = newton_raphson_strategy(
        model_part, time_scheme, linear_solver, convergence_criteria,
        builder_and_solver, solver_settings["max_iterations"].GetInt(),
        solver_settings["compute_reactions"].GetBool(),
        solver_settings["reform_dofs_at_each_step"].GetBool(),
        solver_settings["move_mesh_flag"].GetBool())

    builder_and_solver.SetEchoLevel(solver_settings["echo_level"].GetInt() - 4)
    strategy.SetEchoLevel(solver_settings["echo_level"].GetInt() - 3)
    convergence_criteria.SetEchoLevel(solver_settings["echo_level"].GetInt() - 2)

    if (is_periodic):
        Kratos.Logger.PrintInfo(
            "", "Successfully created periodic solving strategy for " +
            variable.Name() + ".")
    else:
        Kratos.Logger.PrintInfo(
            "", "Successfully created solving strategy for " +
            variable.Name() + ".")

    return strategy


def CreateBuilderAndSolver(linear_solver, is_periodic, communicator):
    if (is_periodic):
        return periodic_block_builder_and_solver(linear_solver, communicator)
    else:
        return block_builder_and_solver(linear_solver, communicator)

def CreateLinearSolver(solver_settings):
    return linear_solver_factory.ConstructSolver(solver_settings)

def CreateScalarScheme(scheme_settings, relaxation_rate, scalar_variable, scalar_variable_rate, relaxed_scalar_variable_rate):
    scheme_type = scheme_settings["scheme_type"].GetString()
    if (scheme_type == "bossak"):
        return bossak_scheme(scheme_settings["scheme_settings"]["alpha_bossak"].GetDouble(),
                             relaxation_rate,
                             scalar_variable,
                             scalar_variable_rate,
                             relaxed_scalar_variable_rate)
    elif (scheme_type == "steady"):
        return steady_scheme(relaxation_rate)
    else:
        raise Exception("Unsupported scheme type. [ scheme_type = \"" + scheme_type + "\" ]")

def InitializePeriodicConditions(model_part, base_model_part, variable):
    properties = model_part.CreateNewProperties(
        model_part.NumberOfProperties() + 1)
    pcu = KratosCFD.PeriodicConditionUtilities(
        model_part, model_part.ProcessInfo[Kratos.DOMAIN_SIZE])
    pcu.AddPeriodicVariable(properties, variable)

    index = model_part.NumberOfConditions()
    for condition in base_model_part.Conditions:
        if condition.Is(Kratos.PERIODIC):
            index += 1
            node_id_list = [node.Id for node in condition.GetNodes()]
            periodic_condition = model_part.CreateNewCondition(
                "PeriodicCondition", index, node_id_list, properties)
            periodic_condition.Set(Kratos.PERIODIC)