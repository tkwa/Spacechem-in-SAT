from ortools.sat.python import cp_model
from game import SpacechemGame, Command, BondDir, Movement
from solution_printer import SolutionPrinter
def run():
    """
    Uses N and H to make NH3.
    """
    N_pattern = [
        [None, None, None, None],
        [None, (7, ""), None, None],
        [None, None, None, None],
        [None, None, None, None],
    ]
    H_pattern = [
        [None, None, None, None],
        [None, (1, ""), None, None],
        [None, None, None, None],
        [None, None, None, None],
    ]
    NH3_pattern = [
        [None, None, None, None],
        [None, (1, "D"), None, None],
        [(1, "R"), (7, "LUR"), (1, "L"), None],
        [None, None, None, None],
    ]
    game = SpacechemGame(T=30, width=10, height=8, n_atom_types=9, max_atoms=4, n_waldos=1, max_bfs=2, n_bonders=4)
    game.check()
    game.make_empty_board_constraint(0)
    game.make_io_constraints(alpha=H_pattern, beta=N_pattern, psi=NH3_pattern)

    game.standard_objective()
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = 3600
    status = solver.Solve(game.model, SolutionPrinter(game, game.width, game.height, print_level='verbose'))

    assert status == cp_model.FEASIBLE or status == cp_model.OPTIMAL, "It Takes Three failed"
    print("Status:", solver.StatusName(status))
    print(f"It Takes Three passed in {solver.WallTime():.3f} s")

if __name__ == "__main__":
    run()