from game import *
from solution_printer import SolutionPrinter

def test_nothing():
    game = SpacechemGame(T=11, width=10, height=8, n_atom_types=2, max_atoms=2, n_waldos=1, max_bfs=5)
    game.check()

    solver = cp_model.CpSolver()
    status = solver.Solve(game.model,None)
    
    assert status == cp_model.FEASIBLE or status == cp_model.OPTIMAL, "Nothing test failed"
    print("Nothing test passed")

def test_atom_BFS():
    game = SpacechemGame(T=5, width=10, height=8, n_atom_types=2, max_atoms=2, n_waldos=1, max_bfs=5)
    game.check()

    game.make_atom_location_constraints(
    {0: [(1, 1, 1),]})
    game.make_waldo_location_constraints(
    [(1, 1, 1),])
    game.model.Add(game.atoms[0].BFS_depth[1] == 0)
    game.model.Add(game.atoms[0].BFS_depth[4] == -1)

    solver = cp_model.CpSolver()
    status = solver.Solve(game.model, None)

    assert status != cp_model.INFEASIBLE, "Atom BFS test failed"
    print("Atom BFS test passed")

def test_atom_grabbed_state():
    game = SpacechemGame(T=11, width=10, height=8, n_atom_types=2, max_atoms=2, n_waldos=1, max_bfs=5)
    game.check()

    game.make_atom_location_constraints(
    {0: [(1, 1, 0),]})
    game.make_waldo_location_constraints(
    [(1, 1, 0),])
    game.model.Add(game.atoms[0].molecule_grabbed[0] == 1)
    game.model.Add(game.atoms[0].molecule_grabbed[8] == 0)

    solver = cp_model.CpSolver()
    status = solver.Solve(game.model, None)

    assert status != cp_model.INFEASIBLE, "Atom molecule_grabbed test failed"
# test_atom_grabbed_state()

def test_waldo_grab_state():
    game = SpacechemGame(T=11, width=10, height=8, n_atom_types=2, max_atoms=2, n_waldos=1, max_bfs=5)
    game.check()

    game.make_atom_location_constraints(
    {0: [(1, 1, 0),]})
    game.make_waldo_location_constraints(
    [(1, 1, 0),])
    game.model.Add(game.waldos[0].grab_active[0] == 1)
    game.model.Add(game.waldos[0].grab_active[8] == 0)

    solver = cp_model.CpSolver()
    status = solver.Solve(game.model, None)

    assert status != cp_model.INFEASIBLE, "Waldo grab_active state change test failed"
    print("Waldo grab_active state change test passed")
# test_waldo_grab_state()

def test_hydrogen_transport():
    game = SpacechemGame(T=9, width=10, height=8, n_atom_types=2, max_atoms=2, n_waldos=1, max_bfs=5)
    game.check()
    game.make_atom_location_constraints(
    {0: [(1, 6, 0), (7, 6, 7)
         ],
     1: [(2, 6, 0), (8, 6, 7),
          ],
     })
    # game.make_waldo_location_constraints(
    # [(1,1, 0)])
    solver = cp_model.CpSolver()
    game.minimize_symbols()
    status = solver.Solve(game.model, SolutionPrinter(game, game.width, game.height))

    assert status != cp_model.INFEASIBLE, "Transport hydrogen test failed"
    assert solver.ObjectiveValue() == 1, f"Transport hydrogen test has wrong number of symbols {solver.ObjectiveValue()}"
    print("Status:", solver.StatusName(status))
    print(f"Time: {solver.WallTime():.3f} s")
    print(f"Stats: {solver.ResponseStats()}")

    print(f"Hydrogen test passed in {solver.WallTime():.3f} s")
# test_hydrogen_transport()

def test_water_input():
    """
    Test that water can be input from input alpha.
    ....
    .HO.
    ..H.
    ....
    """
    game = SpacechemGame(T=3, width=10, height=8, n_atom_types=9, max_atoms=3, n_waldos=1, max_bfs=5)
    game.check()
    game.make_empty_board_constraint(1)
    game.make_input_pattern_constraints(
        pattern=[
        [None, None, None, None],
        [None, (1, "R"), (8, "LD"), None],
        [None, None, (1, "U"), None],
        [None, None, None, None],
        ]
        , input_command=Command.INPUT_ALPHA
    )
    game.make_waldo_location_constraints([(1, 1, 0)])
    game.model.Add(game.atoms[0].x[2] == 1)
    game.model.Add(game.atoms[0].y[2] == 6)
    game.model.Add(game.atoms[0].active[0] == 0)
    game.model.Add(game.atoms[0].active[2] == 1)
    game.model.Add(game.cells[1][6].occupied[2] == 1)
    game.model.Add(game.cells[1][6].occupied[0] == 0)
    game.model.Add(game.cells[1][6].atom_id[2] == 0)
    game.model.Add(game.cells[1][6].atom_at_cell[0][0] == 0)
    game.model.Add(game.cells[1][6].atom_at_cell[2][0] == 1)
    game.model.Add(game.cells[1][6].atom_type[2] == 1)
    game.model.Add(game.cells[2][6].atom_type[2] == 8)
    game.model.Add(game.cells[2][5].atom_type[2] == 1)
    # Should fail test when input alpha cannot be used at t=2
    game.model.Add(game.waldos[0].command[0][Command.INPUT_ALPHA] == 0)

    solver = cp_model.CpSolver()
    game.minimize_symbols()
    status = solver.Solve(game.model, None) # SolutionPrinter(game, game.width, game.height))
    
    assert status != cp_model.INFEASIBLE, "Input water test failed"
    print("Status:", solver.StatusName(status))
    print(f"Water input test passed in {solver.WallTime():.3f} s")

def test_CF4():
    """
    Test that CF4 can be input from input alpha, while other molecules are in the same zone.
    """
    game = SpacechemGame(T=2, width=10, height=8, n_atom_types=9, max_atoms=16, n_waldos=1, max_bfs=5)
    game.check()
    id = 0
    for x in range(0, 4):
        for y in range(4, 8):
            if (x, y) not in [(2, 5), (1, 6), (2, 6), (3, 6), (2, 7)]:
                game.make_atom_location_constraints(
                    {id: [(x, y, 0)]}
                )
                id += 1

    game.make_input_pattern_constraints(
        pattern=[
        [None, None, (9, 'D'), None],
        [None, (9, "R"), (6, "URLD"), (9, 'L')],
        [None, None, (9, "U"), None],
        [None, None, None, None],
        ]
        , input_command=Command.INPUT_ALPHA
    )
    game.model.Add(game.waldos[0].command[1][Command.INPUT_ALPHA] == 1)
    game.model.Add(game.n_active_atoms[1] == 16)

    solver = cp_model.CpSolver()
    game.minimize_symbols()
    status = solver.Solve(game.model, None) # SolutionPrinter(game, game.width, game.height))
    
    assert status != cp_model.INFEASIBLE, "Input CF4 crowded test failed"
    print("Status:", solver.StatusName(status))
    print(f"Input CF4 crowded test passed in {solver.WallTime():.3f} s")

water_pattern =[
        [None, None, None, None],
        [None, (1, "R"), (8, "LD"), None],
        [None, None, (1, "U"), None],
        [None, None, None, None],
        ]

def test_water_transport():
    """
    Test that water can input, then moved from input alpha to output omega.
    ....
    .HO.
    ..H.
    ....
    """
    game = SpacechemGame(T=13, width=10, height=8, n_atom_types=9, max_atoms=3, n_waldos=1, max_bfs=2)
    game.check()
    game.make_empty_board_constraint(0)
    game.make_input_pattern_constraints(
        pattern=[
        [None, None, None, None],
        [None, (1, "R"), (8, "LD"), None],
        [None, None, (1, "U"), None],
        [None, None, None, None],
        ]
        , input_command=Command.INPUT_ALPHA
    )
    game.model.Add(game.cells[8][2].atom_type[12] == 8)
    # game.model.Add(game.cells[8][6].atom_type[8] == 8)
    # game.model.Add(game.waldos[0].command[1][Command.INPUT_ALPHA] == 1)

    solver = cp_model.CpSolver()
    game.minimize_symbols()
    status = solver.Solve(game.model, None) # SolutionPrinter(game, game.width, game.height))
    
    assert status != cp_model.INFEASIBLE, "Transport water test failed"
    print("Status:", solver.StatusName(status))
    print(f"Water transport passed in {solver.WallTime():.3f} s")

def test_input_3_water():
    """Tests that water can be input 3 times."""
    T=9 # minimal time
    game = SpacechemGame(T=T, width=10, height=8, n_atom_types=9, max_atoms=9, n_waldos=1, max_bfs=2)
    game.check()
    game.make_empty_board_constraint(0)
    game.make_input_pattern_constraints(
        pattern=water_pattern
        , input_command=Command.INPUT_ALPHA
    )
    game.make_output_pattern_constraints(pattern=None, output_command=Command.OUTPUT_PSI)
    game.model.Add(sum(game.waldos[0].command[t][Command.INPUT_ALPHA] for t in range(T)) == 3)

    solver = cp_model.CpSolver()
    # game.minimize_symbols()
    status = solver.Solve(game.model, SolutionPrinter(game, game.width, game.height))
    
    assert status != cp_model.INFEASIBLE, "Input 3 water test failed"
    print(f"Status: {solver.StatusName(status)}; {solver.ObjectiveValue()}")
    print(f"Input 3 water passed in {solver.WallTime():.3f} s")

def test_input_output_neon():
    """
    Test that neon can be input, then moved from input alpha to output omega.
    """
    T=10
    game = SpacechemGame(T=T, width=10, height=8, n_atom_types=10, max_atoms=2, n_waldos=1, max_bfs=2)
    game.check()
    game.make_empty_board_constraint(0)
    game.make_atom_location_constraints({0: [(8, 6, 8)]})
    neon_pattern = [
        [None, None, None, None],
        [None, None, (10, ""), None],
        [None, None, None, None],
        [None, None, None, None],
    ]
    game.make_input_pattern_constraints(pattern=neon_pattern, input_command=Command.INPUT_ALPHA)
    game.make_output_pattern_constraints(pattern=neon_pattern, output_command=Command.OUTPUT_PSI)
    game.model.Add(game.waldos[0].command[1][Command.INPUT_ALPHA] == 1)
    game.model.Add(game.waldos[0].command[9][Command.OUTPUT_PSI] == 1)

    solver = cp_model.CpSolver()
    game.minimize_symbols()
    status = solver.Solve(game.model, SolutionPrinter(game, game.width, game.height, print_level='boards'))
    
    assert status != cp_model.INFEASIBLE, "Transport neon test failed"
    print("Status:", solver.StatusName(status))
    print(f"Neon transport passed in {solver.WallTime():.3f} s")

def test_loop_constraint():
    """
    Test that oxygen (O2) can be input and output in a loop.
    """
    oxygen_pattern = [
        [None, None, None, None],
        [None, (8, "R"), (8, "L"), None],
        [None, None, None, None],
        [None, None, None, None],
    ]
    game = SpacechemGame(T=19, width=10, height=8, n_atom_types=9, max_atoms=2, n_waldos=1, max_bfs=2)
    game.check()
    game.make_empty_board_constraint(0)
    game.make_input_pattern_constraints(pattern=oxygen_pattern, input_command=Command.INPUT_ALPHA)
    game.make_output_pattern_constraints(pattern=oxygen_pattern, output_command=Command.OUTPUT_PSI)
    # game.model.Add(game.waldos[0].command[1][Command.INPUT_ALPHA] == 1)
    # Should be doable in a 14-cycle loop
    game.standard_objective()

    solver = cp_model.CpSolver()
    status = solver.Solve(game.model, SolutionPrinter(game, game.width, game.height, print_level='boards'))

    assert status != cp_model.INFEASIBLE, "Loop constraint test failed"
    print("Status:", solver.StatusName(status))
    print(f"Loop constraint passed in {solver.WallTime():.3f} s")

def test_bond_plus():
    """
    Test that O can be bonded using bonders to form O3.
    """
    O_pattern = [
        [None, None, None, None],
        [None, None, (8, ""), None],
        [None, None, None, None],
        [None, None, None, None],
    ]
    O3_pattern = [
        [None, None, None, None],
        [None, (8, "R"), (8, "LD"), None],
        [None, None, (8, "U"), None],
        [None, None, None, None],
    ]
    game = SpacechemGame(T=23, width=10, height=8, n_atom_types=9, max_atoms=3, n_waldos=1, max_bfs=2, n_bonders=4)
    game.check()
    game.make_empty_board_constraint(0)
    game.make_input_pattern_constraints(pattern=O_pattern, input_command=Command.INPUT_ALPHA)
    game.make_output_pattern_constraints(pattern=O3_pattern, output_command=Command.OUTPUT_PSI)
    game.standard_objective()

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = 30
    # game.minimize_symbols()
    status = solver.Solve(game.model, SolutionPrinter(game, game.width, game.height, print_level='boards'))

    assert status == cp_model.FEASIBLE, "Bond+ test failed"
    print("Status:", solver.StatusName(status))
    print(f"Bond+ passed in {solver.WallTime():.3f} s")

def test_bond_minus():
    """
    Test that F2 can be debonded using bonders to form F.
    This is identical to the level "Removing Bonds" in the game.
    The solver's output behavior is different, but the solution found works anyway.
    ┌───────────────────────────────────────┐
    │   │   │   │   │   │   │   │   │   │   │
    │.  │>  │.  │.  │.  │.  │.  │.  │v  │.  │
    │   │ g │ α │W  │   │   │   │ ψ │ d │ g │
    │.  │^  │.  │.  │.  │.  │.  │.  │.  │<  │
    │   │   │   │   │   │   │   │   │ ⊖ │   │
    │.  │.  │.  │.  │.  │.  │.  │.  │.  │.  │
    │   │   │   │   │   │   │   │   │ ψ │   │
    │.  │.  │.  │.  │.  │.  │.  │.  │>  │^  │
    │   │   │   │   │   │   │   │   │   │   │
    │.  │.  │.  │.  │.  │.  │.  │.  │.  │.  │
    │   │   │   │   │   │   │   │   │   │   │
    │.  │.  │.  │.  │.  │.  │.  │.  │.  │.  │
    │   │   │   │   │   │   │   │   │   │   │
    │.  │.  │.  │.  │.  │.  │.  │.  │.  │.  │
    │   │   │   │   │   │   │   │   │   │   │
    │.  │.  │.  │.  │.  │.  │.  │.  │.  │.  │
    └───────────────────────────────────────┘
    """
    F2_pattern = [
        [None, None, None, None],
        [None, (9, "R"), (9, "L"), None],
        [None, None, None, None],
        [None, None, None, None],
    ]
    F_pattern = [
        [None, None, None, None],
        [None, None, (9, ""), None],
        [None, None, None, None],
        [None, None, None, None],
    ]
    game = SpacechemGame(T=23, width=10, height=8, n_atom_types=9, max_atoms=2, n_waldos=1, max_bfs=2, n_bonders=2)
    game.check()
    game.make_empty_board_constraint(0)
    game.make_input_pattern_constraints(pattern=F2_pattern, input_command=Command.INPUT_ALPHA)
    game.make_output_pattern_constraints(pattern=F_pattern, output_command=Command.OUTPUT_PSI)
    # game.model.Add(game.waldos[0].command[1][Command.INPUT_ALPHA] == 1)
    game.standard_objective()

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = 30
    status = solver.Solve(game.model, SolutionPrinter(game, game.width, game.height, print_level='boards'))

    assert status == cp_model.FEASIBLE or status == cp_model.OPTIMAL, "Bond- test failed"
    print("Status:", solver.StatusName(status))
    print(f"Bond- passed in {solver.WallTime():.3f} s")


test_nothing()
test_atom_BFS()
test_water_input()
test_CF4()
test_water_transport()
test_input_3_water()
test_input_output_neon()
test_loop_constraint()

test_bond_plus()
test_bond_minus()


print("All tests passed")