from ortools.sat.python import cp_model
from game import SpacechemGame, Command, BondDir, Movement

class SolutionPrinter(cp_model.CpSolverSolutionCallback):
    """Print intermediate solutions."""

    command_display_dict = {
        Command.NONE: ' ',
        Command.GRAB: 'g',
        Command.DROP: 'd',
        Command.INPUT_ALPHA: 'α',
        Command.INPUT_BETA: 'β',
        Command.OUTPUT_PSI: 'ψ',
        Command.BOND_PLUS: '⊕',
        Command.BOND_MINUS: '⊖',
    }

    def __init__(self, game: SpacechemGame, width, height, print_level='verbose'):
        cp_model.CpSolverSolutionCallback.__init__(self)
        self.model = game.model
        self.__solution_count = 0
        self.game = game
        self.width = width
        self.height = height
        self.print_level = print_level

    def on_solution_callback(self):
        game = self.game
        self.__solution_count += 1
        print(f"Solution {self.__solution_count}")
        for t in range(self.game.T):
            print(f"Time {t}" + (": LOOP START" if t==self.Value(game.t_loop_start) else ": LOOP END" if t==self.Value(game.t_loop_end) else ""))
            # Board state. Each cell is a 4x3 block.
            # Blocks are arranged with bottom left corner at (0, 0)
            # In each block, top left=waldo, top right = command, bottom left=arrow, middle-bottom right=atom type
            state = [[' ' for _ in range(self.width*4)] for _ in range(self.height*3)]
            for waldo in self.game.waldos:
                state[self.Value(waldo.y[t]) * 3][self.Value(waldo.x[t]) * 4] = '@' if self.Value(waldo.grab_active[t]) else 'W'
            for atom in self.game.atoms:
                x, y = self.Value(atom.x[t]), self.Value(atom.y[t])
                if self.Value(atom.active[t]):
                    symbol = SpacechemGame.atom_types_list[self.Value(atom.type[t])]
                    state[self.Value(atom.y[t]) * 3 + 1][self.Value(atom.x[t]) * 4 + 2] = symbol[-1]
                    if len(symbol) == 2:
                        state[self.Value(atom.y[t]) * 3 + 1][self.Value(atom.x[t]) * 4 + 1] = symbol[0]
                    if self.Value(atom.bonds[t][BondDir.L]):
                        state[self.Value(atom.y[t]) * 3 + 1][self.Value(atom.x[t]) * 4 + 1] = '-'
                    if self.Value(atom.bonds[t][BondDir.U]):
                        state[y * 3][x * 4 + 2] = '|'

            for l in self.game.cells:
                for cell in l:
                    for arrow in Movement:
                        if self.Value(cell.arrow[arrow]):
                            state[self.Value(cell.y) * 3 + 1][self.Value(cell.x) * 4] = '>^<v.'[arrow]
                    for command in Command:
                        if self.Value(cell.command[command]):
                            state[self.Value(cell.y) * 3][self.Value(cell.x) * 4 + 1] = self.command_display_dict[command]

            # Print board state encased in double box characters
            print('┌' + '─' * (self.width * 4 - 1) + '┐')
            for y in range(self.height)[::-1]:
                for i in range(2):
                    print('│', end='')
                    print('│'.join([''.join(state[y * 3 + i][x * 4:x * 4 + 3]) for x in range(self.width)]), end='')
                    print('│')
            print('└' + '─' * (self.width * 4 - 1) + '┘')

            # Descriptions of all atoms
            if self.print_level != 'verbose': continue
            for atom in self.game.atoms:
                atom_active = self.Value(atom.active[t])
                if not atom_active:
                    continue
                x, y = self.Value(atom.x[t]), self.Value(atom.y[t])
                cell = self.game.cells[x][y]
                print(f"Atom {atom.id} at {(x,y)}", end=', ')
                if 1 not in [self.Value(atom.movement[t][i]) for i in range(5)]:
                    print(f"Atom {atom.id} has illegal movement!")
                else:
                    print(f"movement {'RULDS'[([self.Value(atom.movement[t][i]) for i in range(5)] + [1]).index(1)]}", end='')
                print(f', bonds {"".join((">^<v"[bond] if self.Value(atom.bonds[t][bond]) else  "") for bond in BondDir)}', end='')
                print(', grabbed' if self.Value(atom.grabbed[t]) else '', end='')
                print(', molecule grabbed' if self.Value(atom.molecule_grabbed[t]) else '', end='')
                print(', BFS depth', self.Value(atom.BFS_depth[t]), end='')
                # print(f', IS at cell' if self.Value(cell.atom_at_cell[t][atom.id]) else ', NOT at cell', end='')
                print(f', cell bonds {"".join((">^<v"[bond] if self.Value(cell.bonds[t][bond]) else  "") for bond in BondDir)}', end='')
                print()
            for waldo in self.game.waldos:
                x, y = self.Value(waldo.x[t]), self.Value(waldo.y[t])
                print(f"Waldo {waldo.id} at {x, y}", end=', ')
                # movement = [self.Value(waldo.movement[t][i]) for i in range(5)]
                # print(movement)
                if 1 not in [self.Value(waldo.movement[t][i]) for i in range(5)]:
                    print(f"Waldo {waldo.id} has illegal movement!")
                else:
                    print(f"movement {'RULDS'[([self.Value(waldo.movement[t][i]) for i in range(5)] + [1]).index(1)]}", end='')
                print(', grabbing' if self.Value(waldo.grab_active[t]) else '', end='')
                commands = [self.Value(waldo.command[t][c]) for c in Command]
                assert sum(commands) >= 1, f"Waldo {waldo.id} has illegal command!"
                print(f', command {Command._member_names_[commands.index(1)]}', end='')
                arrow = [self.Value(waldo.arrow[t][m]) for m in Movement].index(1)
                print(f', arrow {">^<v."[arrow]}')
            n_cells_occupied = sum([self.Value(cell.occupied[t]) for l in self.game.cells for cell in l])
            print(f"{n_cells_occupied} cells occupied by {self.Value(self.game.n_active_atoms[t])} atoms")
            # for l in self.game.cells:
            #     for cell in l:
            #         print(f"Cell {cell.x, cell.y} cell bonds: {''.join(('>^<v'[bond] if self.Value(cell.bonds[t][bond]) else  '') for bond in BondDir)}")
            print()
        print(f"Solution count: {self.solution_count()}")
        print(f"Time: {cp_model.CpSolverSolutionCallback.WallTime(self)}")

    def solution_count(self):
        return self.__solution_count