# %%

from enum import IntEnum
import itertools
from ortools.sat.python import cp_model


class Movement(IntEnum):
    R = 0
    U = 1
    L = 2
    D = 3
    # ROTATE_CCW = 4
    # ROTATE_CW = 5
    STALL = 4

class BondDir(IntEnum):
    R = 0
    U = 1
    L = 2
    D = 3


class Atom():
    """
    Represents an atom id. Atom ids can be reused.
    T is max timesteps.
    """
    def __init__(self, model:cp_model.CpModel, T:int, id:int, reactor_height, reactor_width, n_atom_types, max_bfs):
        rT = range(T)
        self.model = model
        self.T = T
        self.id = id
        self.n_atom_types = n_atom_types

        self.active = [model.NewBoolVar(f'atom_{id}_active_{t}') for t in rT]
        self.type = [model.NewIntVar(0, n_atom_types, f'atom{id}_type_{t}') for t in rT]
        # x and y are stored as integers
        self.x = [model.NewIntVar(0, reactor_width -1, f'atom_{id}_x_{t}') for t in rT]
        self.y = [model.NewIntVar(0, reactor_height-1, f'atom_{id}_y_{t}') for t in rT]
        # self.xy = [iv(0, width * height - 1, f'atom_{id}_xy_{t}') for t in rT]
        # boolean version of x and y
        self.xb = [[model.NewBoolVar(f'atom_{id}_x_{x}_{t}') for x in range(reactor_width)] for t in rT]
        self.yb = [[model.NewBoolVar(f'atom_{id}_y_{y}_{t}') for y in range(reactor_height)] for t in rT]
        self.movement = [[model.NewBoolVar(f'atom_{id}_move_{movement}_{t}') for movement in Movement] for t in rT]
        # self.molecule_id = [None]*T
        self.bonds = [[model.NewBoolVar(f'atom_{id}_bond_{bond_dir}_{t}') for bond_dir in BondDir] for t in rT]
        # This will be extended to whether the atom is grabbed by each waldo
        self.grabbed = [model.NewBoolVar(f'atom_{id}_grab_{t}') for t in rT]
        # Distance from the waldo, if grabbed. Otherwise -1.
        self.BFS_depth = [model.NewIntVar(-1, max_bfs, f'atom_{id}_BFS_{t}') for t in rT]
        self.molecule_grabbed = [model.NewBoolVar(f'atom_{id}_molecule_grab_{t}') for t in rT]

    def check(self, t:int):
        model = self.model
        # Check integer position vs boolean position using Element
        # model.AddElement(self.x[t], self.xb[t], 1)
        # model.AddElement(self.y[t], self.yb[t], 1)
        # model.AddExactlyOne(self.xb[t])
        # model.AddExactlyOne(self.yb[t])
        # xy position for AddElement interface with cells
        # model.Add(self.xy[t] == self.x[t] + self.y[t] * self.width).OnlyEnforceIf(self.active[t])


        # Only grabbed atoms can move
        model.Add(self.movement[t][Movement.STALL] == 1).OnlyEnforceIf(self.molecule_grabbed[t].Not())
        # Atom movement.
        model.AddExactlyOne(self.movement[t])
        if t < self.T-1:
            model.Add(self.x[t+1] == self.x[t] + 1).OnlyEnforceIf(self.movement[t][Movement.R])
            model.Add(self.x[t+1] == self.x[t] - 1).OnlyEnforceIf(self.movement[t][Movement.L])
            model.Add(self.y[t+1] == self.y[t] + 1).OnlyEnforceIf(self.movement[t][Movement.U])
            model.Add(self.y[t+1] == self.y[t] - 1).OnlyEnforceIf(self.movement[t][Movement.D])
            for movement in Movement.R, Movement.L, Movement.STALL:
                model.Add(self.y[t+1] == self.y[t]).OnlyEnforceIf(self.movement[t][movement])
            for movement in Movement.U, Movement.D, Movement.STALL:
                model.Add(self.x[t+1] == self.x[t]).OnlyEnforceIf(self.movement[t][movement])

            # Active atoms cannot change type.
            model.Add(self.type[t+1] == self.type[t]).OnlyEnforceIf(self.active[t])

            # Atom bonds are preserved (except bond command which come later). This is cheaper than checking cells.
            for bond_dir in BondDir:
                model.Add(self.bonds[t+1][bond_dir] == self.bonds[t][bond_dir]).OnlyEnforceIf(self.active[t])

        model.Add(self.BFS_depth[t] == -1).OnlyEnforceIf(self.molecule_grabbed[t].Not())
        model.Add(self.BFS_depth[t] >= 0).OnlyEnforceIf(self.molecule_grabbed[t]) # this is the problem...

# %%

class Command(IntEnum):
    NONE = 0
    GRAB = 1
    DROP = 2
    # GRABDROP = 3
    # ROTATE_CCW = 4
    # ROTATE_CW = 5
    # BOND_PLUS = 6
    # BOND_MINUS = 7
    INPUT_ALPHA = 3
    # INPUT_BETA = 9
    # OUTPUT_PSI = 4
    # OUTPUT_OMEGA = 11
    # Still need to implement: flip flops, sense, fuse/split

class Cell():
    def __init__(self, model:cp_model.CpModel, T, height, width, x:int, y:int, n_atom_types, max_atoms, max_bfs):
        self.model = model
        self.height = height
        self.width = width
        self.x = x
        self.y = y
        # Whether the cell has an active atom with matching x,y.
        self.occupied = [model.NewBoolVar(f'cell_{x}_{y}_occupied_{t}') for t in range(T)]
        self.atom_at_cell = [[model.NewBoolVar(f'cell_{x}_{y}_atom_{id} at cell {t}') for id in range(max_atoms)] for t in range(T)]
        # atom_id is undefined if no atom is present; atom_type is 0 if no atom is present.
        self.atom_id = [model.NewIntVar(0, max_atoms-1, f'cell_{x}_{y}_atom_id_{t}') for t in range(T)]
        self.atom_type = [model.NewIntVar(0, n_atom_types, f'cell_{x}_{y}_atom_type_{t}') for t in range(T)]
        self.waldo_at_cell = [model.NewBoolVar(f'cell_{x}_{y}_waldo_{t}') for t in range(T)]
        self.waldo_grabbing = [model.NewBoolVar(f'cell_{x}_{y}_waldo_grabbing_{t}') for t in range(T)]
        self.command = [model.NewBoolVar(f'cell_{x}_{y}_command_{c}') for c in Command]
        self.arrow = [model.NewBoolVar(f'cell_{x}_{y}_arrow_{v}') for v in Movement]
        self.bonds = [[model.NewBoolVar(f'cell_{x}_{y}_bond_{bond_dir}_{t}') for bond_dir in BondDir] for t in range(T)]
        self.BFS_depth = [model.NewIntVar(-1, max_bfs, f'cell_{x}_{y}_bfs_depth_{t}') for t in range(T)]
        self.BFS_depth_ge_1 = [model.NewBoolVar(f'cell_{x}_{y}_indirectly grabbed_{t}') for t in range(T)]
        self.BFS_parent_dirs = [[model.NewBoolVar(f'cell_{x}_{y}_parent_direction_{bond_dir}_{t}') for bond_dir in BondDir] for t in range(T)]

    def check(self, t:int):
        self.model.AddImplication(self.waldo_grabbing[t], self.occupied[t])
        self.model.AddImplication(self.waldo_grabbing[t], self.waldo_at_cell[t])

        # At most one arrow (or None/Stall) and exactly one command
        self.model.AddExactlyOne(self.arrow)
        self.model.AddExactlyOne(self.command)

        # Atoms
        self.model.AddExactlyOne(self.atom_at_cell[t] + [self.occupied[t].Not()])
        self.model.Add(self.atom_type[t] == 0).OnlyEnforceIf(self.occupied[t].Not())
        self.model.Add(self.atom_type[t] != 0).OnlyEnforceIf(self.occupied[t])

        # Bonding
        # A cell can only have bonds if it is occupied
        for bond_dir in BondDir:
            self.model.Add(self.bonds[t][bond_dir] == 0).OnlyEnforceIf(self.occupied[t].Not())

        # BFS depth
        # any square not equal to 1 or 0 must have at least one neighbor with value v-1
        # Unknown depth if occupied but not waldo
        self.model.Add(self.BFS_depth[t] >= 1).OnlyEnforceIf(self.BFS_depth_ge_1[t])
        self.model.Add(self.BFS_depth[t] <  1).OnlyEnforceIf(self.BFS_depth_ge_1[t].Not())
        self.model.Add(self.BFS_depth[t] == 0).OnlyEnforceIf(self.waldo_at_cell[t], self.waldo_grabbing[t])
        self.model.Add(self.BFS_depth[t] != 0).OnlyEnforceIf(self.waldo_at_cell[t].Not())
        self.model.Add(self.BFS_depth[t] != 0).OnlyEnforceIf(self.waldo_grabbing[t].Not())
        self.model.Add(self.BFS_depth[t] == -1).OnlyEnforceIf(self.occupied[t].Not())
        possible_parents = [] # to check if any of these are equal to BFS_depth[t]-1
        # If molecule has BFS depth k>0, it must have a parent that it's bonded to with depth k-1
        # Only check directions that have bonds
        for bond_dir in BondDir:
            self.model.AddImplication(self.BFS_parent_dirs[t][bond_dir], self.bonds[t][bond_dir])
        if self.x > 0:
            possible_parents.append(self.BFS_parent_dirs[t][0])
        if self.x < self.width - 1:
            possible_parents.append(self.BFS_parent_dirs[t][1])
        if self.y > 0:
            possible_parents.append(self.BFS_parent_dirs[t][2])
        if self.y < self.height - 1:
            possible_parents.append(self.BFS_parent_dirs[t][3])
        self.model.AddBoolOr(possible_parents).OnlyEnforceIf(self.BFS_depth_ge_1)
        self.model.Add(self.BFS_depth[t] <= 0).OnlyEnforceIf([var.Not() for var in possible_parents])

# %%
class Waldo():
    def __init__(self, model:cp_model.CpModel, T, id:int, height, width, max_atoms):
        bv = model.NewBoolVar
        iv = model.NewIntVar
        self.T = T
        self.model = model
        self.id = id
        self.max_atoms = max_atoms
        self.x = [iv(0, width -1, f'waldo_{id}_x_{t}') for t in range(T)]
        self.y = [iv(0, height-1, f'waldo_{id}_y_{t}') for t in range(T)]
        self.movement = [[bv(f'waldo_{id}_move_{movement}_{t}') for movement in Movement] for t in range(T)]
        self.grab_active = [bv(f'waldo_{id}_grabbing_{t}') for t in range(T)]
        self.grabbed_atom = [[bv(f'waldo_{id}_grabbed_atom{atom_id}_{t}') for atom_id in range(max_atoms)] for t in range(T)]
        self.command = [[bv(f'waldo_{id}_command_{command}_{t}') for command in Command] for t in range(T)]
        self.arrow = [[bv(f'waldo_{id}_arrow_{m}_{t}') for m in Movement] for t in range(T)]

    def check(self, t:int):
        model = self.model
        model.AddExactlyOne(self.movement[t])
        model.AddExactlyOne(self.command[t])
        # In SpaceChem, the start command takes up cycle 0.
        model.Add(self.command[0][Command.NONE] == 1)
        model.AddExactlyOne(self.arrow[t])
        # Movement
        if t < self.T - 1:
            model.Add(self.x[t+1] == self.x[t] + 1).OnlyEnforceIf(self.movement[t][Movement.R])
            model.Add(self.x[t+1] == self.x[t] - 1).OnlyEnforceIf(self.movement[t][Movement.L])
            model.Add(self.y[t+1] == self.y[t] + 1).OnlyEnforceIf(self.movement[t][Movement.U])
            model.Add(self.y[t+1] == self.y[t] - 1).OnlyEnforceIf(self.movement[t][Movement.D])
            for movement in Movement.R, Movement.L, Movement.STALL:
                model.Add(self.y[t+1] == self.y[t]).OnlyEnforceIf(self.movement[t][movement])
            for movement in Movement.U, Movement.D, Movement.STALL:
                model.Add(self.x[t+1] == self.x[t]).OnlyEnforceIf(self.movement[t][movement])
        # Movement matches arrow
        if t > 0:
            for movement in Movement:
                model.Add(self.movement[t][movement] == self.movement[t-1][movement]).OnlyEnforceIf(self.arrow[t][Movement.STALL])
        for movement in Movement.U, Movement.D, Movement.L, Movement.R:
            model.Add(self.movement[t][movement] == 1).OnlyEnforceIf(self.arrow[t][movement])

        # Grabbing.
        # Case where waldo grabs an atom is handled in SpacechemGame.check().
        # We can't use AddExactlyOne with OnlyEnforceIf, so we use an equivalent expression.
        model.Add(self.grab_active[0] == 0)
        model.AddExactlyOne(self.grabbed_atom[t] + [self.grab_active[t].Not()])
        model.Add(self.grab_active[t] == 0).OnlyEnforceIf(self.command[t][Command.DROP])
        if t > 0:
            pass
            model.Add(self.grab_active[t] == self.grab_active[t-1]).OnlyEnforceIf(self.command[t][Command.NONE])
            model.Add(self.command[t][Command.DROP] == 1).OnlyEnforceIf(self.grab_active[t-1], self.grab_active[t].Not())
            model.Add(self.command[t][Command.GRAB] == 1).OnlyEnforceIf(self.grab_active[t-1].Not(), self.grab_active[t])



        
            
# %%


class SpacechemGame():
    """
    Variables are an object's state at time t, BEFORE atoms move to time t+1.
    """
    atom_types_list = ['?', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne']

    def __init__(self, T:int, width:int=10, height:int=8, n_atom_types:int=100, max_atoms:int=80, n_waldos:int=1, max_bfs:int=10):
        """
        max_bfs is the maximum distance from an atom to a grabbing waldo.
        """
        self.model = cp_model.CpModel()
        self.T = T
        self.width = width
        self.height = height
        self.n_atom_types = n_atom_types
        self.max_atoms = max_atoms
        self.atoms = [Atom(self.model, T, id, height, width, n_atom_types=n_atom_types, max_bfs=max_bfs) for id in range(max_atoms)]
        # Number of active atoms. Atoms become active on the cycle they are input, and inactive on the cycle they are output.
        self.n_active_atoms = [self.model.NewIntVar(0, max_atoms, f'n_active_atoms_{t}') for t in range(T)]
        self.cells:list[list[Cell]] = [[Cell(self.model, T, height, width, x, y, n_atom_types, max_atoms, max_bfs) for y in range(height)] for x in range(width)]
        self.waldos = [Waldo(self.model, T, id, height, width, max_atoms) for id in range(n_waldos)]
        # self.bonders = [[None for _ in range(height)] for _ in range(width)]
    
        assert n_waldos == 1, "Only one waldo is supported for now"

    def atom(self, id):
        return self.atoms[id]
    
    def minimize_symbols(self):
        # Minimizes the non-NONE commands and arrows across all cells
        commands = [cell.command[Command.NONE].Not() for l in self.cells for cell in l]
        arrows = [cell.arrow[Movement.STALL].Not() for l in self.cells for cell in l]
        self.model.Minimize(sum(commands + arrows))
    
    def check(self):
        m = self.model
        for t in range(self.T):
            # Check atoms, cells, and waldos
            for atom in self.atoms:
                atom.check(t)
            for x in range(self.width):
                for y in range(self.height):
                    self.cells[x][y].check(t)
            for waldo in self.waldos:
                waldo.check(t)
            # Atoms cannot become active or inactive (except input alpha and output psi)
            if t > 0:
                for atom in self.atoms:
                    m.Add(atom.active[t] <= atom.active[t-1]).OnlyEnforceIf(self.waldos[0].command[t][Command.INPUT_ALPHA].Not())
                    m.Add(atom.active[t] >= atom.active[t-1]) # .OnlyEnforceIf(self.waldos[0].command[t][Command.OUTPUT_PSI].Not())
                self.model.Add(self.n_active_atoms[t] == self.n_active_atoms[t-1]).OnlyEnforceIf(self.waldos[0].command[t][Command.INPUT_ALPHA].Not())
            # Whether an atom is active is determined by its ID.
            for id, atom in enumerate(self.atoms):
                m.Add(atom.id < self.n_active_atoms[t]).OnlyEnforceIf(atom.active[t])
                m.Add(atom.id >= self.n_active_atoms[t]).OnlyEnforceIf(atom.active[t].Not())

            # No two active atoms have the same position
            for i, atom1 in enumerate(self.atoms):
                for atom2 in self.atoms[i+1:]:
                    m.Add(atom1.y[t] * self.width + atom1.x[t] != atom2.y[t] * self.width + atom2.x[t]).OnlyEnforceIf([atom1.active[t], atom2.active[t]])
            
            
            # Check atoms against waldos
            for id, atom in enumerate(self.atoms):
                assert atom.id == id
                for waldo in self.waldos:
                    for atom in self.atoms:
                        waldo_grabbing_atom = waldo.grabbed_atom[t][atom.id]
                        waldo_grabbing_molecule = atom.molecule_grabbed[t]
                        m.AddImplication(waldo_grabbing_atom, waldo_grabbing_molecule)

                        m.AddImplication(waldo_grabbing_molecule, atom.active[t])
                        m.Add(atom.x[t] == waldo.x[t]).OnlyEnforceIf(waldo_grabbing_atom)
                        m.Add(atom.y[t] == waldo.y[t]).OnlyEnforceIf(waldo_grabbing_atom)
                        #TODO a waldo can grab an atom moving any direction
                        for movement in Movement:
                            m.Add(atom.movement[t][movement] == waldo.movement[t][movement]).OnlyEnforceIf(waldo_grabbing_molecule)

                        m.AddImplication(atom.grabbed[t], waldo_grabbing_atom)
                        m.AddImplication(waldo_grabbing_atom, atom.grabbed[t])

            # AddElement method: 5.3 seconds
            # for id, atom in enumerate(self.atoms):
            #     assert atom.id == id
            #     # remember, atom xy is only constrained if atom is active
            #     m.AddElement(atom.xy[t], [cell.atom_id[t] for cell in self.cells_flat], atom.id)
            #     m.AddElement(atom.xy[t], [cell.atom_type[t] for cell in self.cells_flat], atom.type[t])
            #     m.AddElement(atom.xy[t], [cell.occupied[t] for cell in self.cells_flat], 1)
            # for l in self.cells:
            #     for cell in l:
            #         # If a cell is unoccupied, its atom id can be assigned to the dummy atom
            #         # what's keeping two atoms from being in the same cell? this must be a separate constraint
            #         dummy_x = m.NewIntVar(cell.x, cell.x, f'dummy_x_{cell.x}_{cell.y}')
            #         dummy_y = m.NewIntVar(cell.y, cell.y, f'dummy_y_{cell.x}_{cell.y}')
            #         m.AddElement(cell.atom_id[t], [atom.x[t] for atom in self.atoms] + [dummy_x], cell.x)
            #         m.AddElement(cell.atom_id[t], [atom.y[t] for atom in self.atoms] + [dummy_y], cell.y)
            # Check bonds between cells
            for x in range(self.width):
                for y in range(self.height):
                    cell:Cell = self.cells[x][y]
                    m.Add(cell.bonds[t][BondDir.L] == (self.cells[x-1][y].bonds[t][BondDir.R] if x > 0 else 0))
                    m.Add(cell.bonds[t][BondDir.R] == (self.cells[x+1][y].bonds[t][BondDir.L] if x < self.width - 1 else 0))
                    m.Add(cell.bonds[t][BondDir.D] == (self.cells[x][y-1].bonds[t][BondDir.U] if y > 0 else 0))
                    m.Add(cell.bonds[t][BondDir.U] == (self.cells[x][y+1].bonds[t][BondDir.D] if y < self.height - 1 else 0))

            # Check BFS depth between cells
            for x in range(self.width):
                for y in range(self.height):
                    cell:Cell = self.cells[x][y]
                    # check parents have BFS depth k-1 if cell has BFS depth k
                    if x > 0:
                        m.Add(cell.BFS_depth[t] == self.cells[x-1][y].BFS_depth[t] + 1).OnlyEnforceIf(cell.BFS_parent_dirs[t][BondDir.L])
                    if x < self.width - 1:
                        m.Add(cell.BFS_depth[t] == self.cells[x+1][y].BFS_depth[t] + 1).OnlyEnforceIf(cell.BFS_parent_dirs[t][BondDir.R])
                    if y > 0:
                        m.Add(cell.BFS_depth[t] == self.cells[x][y-1].BFS_depth[t] + 1).OnlyEnforceIf(cell.BFS_parent_dirs[t][BondDir.D])
                    if y < self.height - 1:
                        m.Add(cell.BFS_depth[t] == self.cells[x][y+1].BFS_depth[t] + 1).OnlyEnforceIf(cell.BFS_parent_dirs[t][BondDir.U])

            # Check atoms against cells (expensive; uses 2/3 of the time!)
            # Atom-at-cell method: 3.3 seconds
            for l in self.cells:
                for cell in l:
                    for id, atom in enumerate(self.atoms):
                        assert atom.id == id
                        atom_at_cell_x = m.NewBoolVar(f'atom_{id}_at_cell_{cell.x}_{cell.y}_x_{t}')
                        atom_at_cell_y = m.NewBoolVar(f'atom_{id}_at_cell_{cell.x}_{cell.y}_y_{t}')
                        atom_at_cell = cell.atom_at_cell[t][id]
                        m.Add(atom.x[t] == cell.x).OnlyEnforceIf(atom_at_cell_x)
                        m.Add(atom.x[t] != cell.x).OnlyEnforceIf(atom_at_cell_x.Not())
                        m.Add(atom.y[t] == cell.y).OnlyEnforceIf(atom_at_cell_y)
                        m.Add(atom.y[t] != cell.y).OnlyEnforceIf(atom_at_cell_y.Not())
                        m.AddBoolAnd([atom_at_cell_x, atom_at_cell_y, atom.active[t]]).OnlyEnforceIf(atom_at_cell)
                        m.Add(atom_at_cell == 1).OnlyEnforceIf(atom_at_cell_x, atom_at_cell_y, atom.active[t])
                        # All this is needed for the cell.atom_id, cell.atom_type, and cell.occupied variables, plus checking bonds
                        m.Add(cell.atom_id[t] == atom.id).OnlyEnforceIf(atom_at_cell)
                        m.Add(cell.atom_type[t] == atom.type[t]).OnlyEnforceIf(atom_at_cell)
                        # Check bonds. Bonds can potentially be optimized by a factor of 2.
                        for bond in BondDir:
                            m.AddImplication(atom.bonds[t][bond], cell.bonds[t][bond]).OnlyEnforceIf(atom_at_cell)
                            m.AddImplication(cell.bonds[t][bond], atom.bonds[t][bond]).OnlyEnforceIf(atom_at_cell)
                        # If there is no atom at the cell, BFS depth is undefined
                        m.Add(cell.BFS_depth[t] == atom.BFS_depth[t]).OnlyEnforceIf(atom_at_cell)

            # Check waldos against cells (for commands)
            waldo = self.waldos[0]
            for l in self.cells:
                for cell in l:
                    waldo_at_cell_x = m.NewBoolVar(f'waldo_at_cell_{cell.x}_{cell.y}_x_{t}')
                    waldo_at_cell_y = m.NewBoolVar(f'waldo_at_cell_{cell.x}_{cell.y}_y_{t}')
                    waldo_at_cell = cell.waldo_at_cell[t]
                    m.Add(waldo.x[t] == cell.x).OnlyEnforceIf(waldo_at_cell_x)
                    m.Add(waldo.x[t] != cell.x).OnlyEnforceIf(waldo_at_cell_x.Not())
                    m.Add(waldo.y[t] == cell.y).OnlyEnforceIf(waldo_at_cell_y)
                    m.Add(waldo.y[t] != cell.y).OnlyEnforceIf(waldo_at_cell_y.Not())
                    m.AddBoolAnd([waldo_at_cell_x, waldo_at_cell_y]).OnlyEnforceIf(waldo_at_cell)
                    m.AddBoolOr([waldo_at_cell_x.Not(), waldo_at_cell_y.Not()]).OnlyEnforceIf(waldo_at_cell.Not())
                    # Match waldo command to cell command
                    for command in Command:
                        m.Add(waldo.command[t][command] == 1).OnlyEnforceIf(waldo_at_cell, cell.command[command])
                    for arrow in Movement:
                        m.Add(waldo.arrow[t][arrow] == 1).OnlyEnforceIf(waldo_at_cell, cell.arrow[arrow])

                    # Cell says waldo grabbing iff waldo is grabbing
                    m.Add(waldo.grab_active[t] == 1).OnlyEnforceIf(waldo_at_cell, cell.waldo_grabbing[t])
                    m.Add(waldo.grab_active[t] == 0).OnlyEnforceIf(waldo_at_cell, cell.waldo_grabbing[t].Not())
                    # Waldo GRAB command grabs an atom iff it's at the cell
                    m.Add(waldo.grab_active[t] == cell.occupied[t]).OnlyEnforceIf(waldo_at_cell, waldo.command[t][Command.GRAB])


        

    def run(self, waldo_positions):
        m=self.model
        t = 0
        # TODO set waldo positions

        for t in range(1, self.T):
            # waldo actions
            # for waldo in self.waldos:
            #     x, y = waldo.x[t], waldo.y[t]
            #     waldo.command[t] = self.cells[x, y].command
            #     # waldo actions
            #     match waldo.command[t]:
            #         case Command.NONE:
                        
            #             self.model.Add()
            #             pass # waldos continue moving
            #         case Command.GRAB:
            #             # Constraint: while waldo grabbing, any atom at waldo position is grabbed
            #             # Grabs persist until a drop.
            #             waldo.grabbing[t] = True
            #         case Command.DROP:
            #             waldo.grabbing[t] = False
            #         case Command.GRABDROP:
            #             waldo.grabbing[t] = not waldo.grabbing[t-1]
            #         case Command.ROTATE_CCW:
            #             # TODO rotate atoms
            #             pass
            #         case Command.ROTATE_CW:
            #             pass
            #         case Command.BOND_PLUS:
            #             # Bonds atoms iff they are on bonders, adjacent, and not at max bonds
            #             # TODO bond atoms
            #             pass
            #         case Command.BOND_MINUS:
            #             # Unbonds atoms iff they are on bonders, adjacent, and bonded
            #             # TODO un-bond atoms
            #             pass
            #         case Command.INPUT_ALPHA:
            #             # Makes an input alpha appear
            #             pass
            #         case Command.INPUT_BETA:
            #             pass
            #         case Command.OUTPUT_PSI:
            #             pass
            #         case Command.OUTPUT_OMEGA:
            #             pass
            #         case _:
            #             pass
            #             #raise Exception("Invalid waldo command")
                    
            #     # move waldos
            #     match waldo.movement[t]:
            #         case Waldo.Movement.R:
            #             waldo.x[t+1] = waldo.x[t] + 1
            #             waldo.y[t+1] = waldo.y[t]
            #         case Waldo.Movement.U:
            #             waldo.x[t+1] = waldo.x[t]
            #             waldo.y[t+1] = waldo.y[t] + 1
            #         case Waldo.Movement.L:
            #             waldo.x[t+1] = waldo.x[t] - 1
            #             waldo.y[t+1] = waldo.y[t]
            #         case Waldo.Movement.D:
            #             waldo.x[t+1] = waldo.x[t]
            #             waldo.y[t+1] = waldo.y[t] - 1
            #         case Waldo.Movement.STALL:
            #             waldo.x[t+1] = waldo.x[t]
            #             waldo.y[t+1] = waldo.y[t]
            #         case _:
            #             pass
            #             # raise Exception("Invalid waldo movement")


            # move atoms
            # Checking done in `check` method of atoms
            for atom in self.atoms:
                if not atom.active:
                    continue
                if atom.movement is not Atom.Movement.NONE:
                    assert atom.molecule_grabbed[t]
                match atom.movement[t]:
                    case Atom.Movement.R:
                        atom.x[t+1] = atom.x[t] + 1
                        atom.y[t+1] = atom.y[t]
                    case Atom.Movement.U:
                        atom.x[t+1] = atom.x[t]
                        atom.y[t+1] = atom.y[t] + 1
                    case Atom.Movement.L:
                        atom.x[t+1] = atom.x[t] - 1
                        atom.y[t+1] = atom.y[t]
                    case Atom.Movement.D:
                        atom.x[t+1] = atom.x[t]
                        atom.y[t+1] = atom.y[t] - 1
                    # case Atom.Movement.ROTATE_CCW:
                    #     atom.rotate_center_x[t] = self.molecules[atom.molecule_id].center_x[t]
                    #     atom.rotate_center_y[t] = self.molecules[atom.molecule_id].center_y[t]
                    #     pass
                    # case Atom.Movement.ROTATE_CW:
                    #     pass
                    case Atom.Movement.NONE:
                        pass
                    case _:
                        pass
                        # raise Exception("Invalid atom movement")

            # Collision checking
            # TODO this is O(atoms^2), we can do O(squares)
            for atom1, atom2 in itertools.combinations(self.atoms, 2):
                if not atom1.active or not atom2.active:
                    continue
                assert not(atom1.x[t] == atom2.x[t] and atom1.y[t] == atom2.y[t])


    def make_atom_location_constraints(self, atoms_dict: dict[int, list[tuple[int, int, int]]]):
        """
        Takes a dictionary of atom ids to a list of (x, y, t) tuples.
        """
        for atom_id, tuples in atoms_dict.items():
            for x, y, t in tuples:
                self.model.Add(self.atoms[atom_id].x[t] == x)
                self.model.Add(self.atoms[atom_id].y[t] == y)
                self.model.Add(self.atoms[atom_id].active[t] == True)

    def make_empty_board_constraint(self, t:int=0):
        for x in range(self.width):
            for y in range(self.height):
                self.model.Add(self.cells[x][y].occupied[t] == False)

    def draw_atom_positions(self, board_state: str, t:int=0):
        """
        Takes a string representation of the board state, and enforces those constraints.
        """
        lines = board_state.splitlines()
        assert len(lines) == self.height
        for y, line in enumerate(lines)[::-1]:
            assert len(line) == self.width * 2
            for x in len(self.width):
                sym = line[x*2:x*2+2].strip()
                if sym:
                    atom_type = SpacechemGame.atom_types_list.index(sym)
                    self.model.Add(self.cells[x][y].atom_type[t] == atom_type)

    def make_input_pattern_constraints(self, pattern:list, input_command: Command.INPUT_ALPHA):
        """
        Takes a 4x4 list of atom types, and generates constraints that enforce that pattern.
        At any time t, input can be activated. If so,
        - n new atoms must be made active at time t
        - for each filled cell in the pattern, the cell must have the corresponding atom type at time t
        - the cells are filled with new atom ids

        Pattern format: a 4x4 list where each element is (atom_type, (bondR, bondU, bondL, bondD)) or None
        Increasing y is up
        """
        pattern = pattern[::-1]
        n_new_atoms = sum(1 for row in pattern for atom_spec in row if atom_spec is not None)
        x_offset = 0
        y_offset = 4 if input_command == Command.INPUT_ALPHA else 0
        atoms_dict = {(x + x_offset, y + y_offset): atom for y, row in enumerate(pattern) for x, atom in enumerate(row) if atom is not None}
        assert self.max_atoms >= n_new_atoms

        for t in range(1, self.T):
            input_at_t:cp_model.IntVar = self.waldos[0].command[t][input_command]

            self.model.Add(self.n_active_atoms[t] == self.n_active_atoms[t-1] + n_new_atoms).OnlyEnforceIf(input_at_t)

            for x, y in atoms_dict:
                atom_type, bond_string = atoms_dict[(x,y)]
                bonds = [bond_chr in bond_string for bond_chr in "RULD"]
                cell = self.cells[x][y]
                self.model.Add(cell.occupied[t] == True).OnlyEnforceIf(input_at_t)
                self.model.Add(cell.atom_type[t] == atom_type).OnlyEnforceIf(input_at_t)
                self.model.Add(self.n_active_atoms[t-1] <= cell.atom_id[t]).OnlyEnforceIf(input_at_t)
                self.model.Add(cell.atom_id[t] < self.n_active_atoms[t]).OnlyEnforceIf(input_at_t)
                for bond_dir in BondDir:
                    self.model.Add(cell.bonds[t][bond_dir] == bonds[bond_dir]).OnlyEnforceIf(input_at_t)
        


    def make_waldo_location_constraints(self, waldo_locations: list[tuple[int, int, int]]):
        """
        Takes a list of (x, y, t) tuples.
        """
        for x, y, t in waldo_locations:
            self.model.Add(self.waldos[0].x[t] == x)
            self.model.Add(self.waldos[0].y[t] == y)



# %%

class SolutionPrinter(cp_model.CpSolverSolutionCallback):
    """Print intermediate solutions."""

    command_display_dict = {
        Command.NONE: ' ',
        Command.GRAB: 'g',
        Command.DROP: 'd',
        Command.INPUT_ALPHA: 'α'
        # Command.BOND_PLUS: '+',
    }

    def __init__(self, game: SpacechemGame, width, height):
        cp_model.CpSolverSolutionCallback.__init__(self)
        self.model = game.model
        self.__solution_count = 0
        self.game = game
        self.width = width
        self.height = height

    def on_solution_callback(self):
        game = self.game
        self.__solution_count += 1
        print(f"Solution {self.__solution_count}")
        for t in range(self.game.T):
            print(f"Time {t}")
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
            print(f"{n_cells_occupied} cells occupied")
            # for l in self.game.cells:
            #     for cell in l:
            #         print(f"Cell {cell.x, cell.y} cell bonds: {''.join(('>^<v'[bond] if self.Value(cell.bonds[t][bond]) else  '') for bond in BondDir)}")
            print()

    def solution_count(self):
        return self.__solution_count


# %%

if __name__ == '__main__':
    game = SpacechemGame(T=9, width=10, height=8, n_atom_types=100, max_atoms=4, n_waldos=1, max_bfs=5)
    game.check()

    print(game.model.Validate())

    # game.model.Maximize(sum([game.atoms[id].active[t] for id in range(len(game.atoms)) for t in range(game.T)]))
    game.model.Add(game.atoms[0].type[8] == 1)
    game.model.Add(game.atoms[1].type[8] == 1)

    game.make_atom_location_constraints(
        {0: [(1, 1, 0), (3,3, -1)
            ],
        # 1: [(1, 2, 0), (3, 2, -1)],
        }
    )

    game.make_waldo_location_constraints([
        (1, 1, 0),
    ])

    # game.model.Add(game.cells[1][1].waldo_at_cell[1] == 0)

    game.model.Add(game.waldos[0].x[0] == game.waldos[0].x[-1])
    game.model.Add(game.waldos[0].y[0] == game.waldos[0].y[-1])
    # game.model.Add(game.waldos[0].movement[0][Movement.R] == 1)
    # game.model.Add(game.waldos[0].movement[-1][Movement.R] == 1)
    game.model.Add(game.waldos[0].grab_active[0] == 0)
    game.model.Add(game.waldos[0].grab_active[-1] == 0)

    # game.model.Add(game.atoms[2].active[1] == 0)
    game.minimize_symbols()
    print(f"Done building model")
    # Solve model
    solver = cp_model.CpSolver()
    status = solver.Solve(game.model, SolutionPrinter(game, game.width, game.height))
    print("Status:", solver.StatusName(status))
    print(f"Time: {solver.WallTime():.3f} s")
    print(f"Stats: {solver.ResponseStats()}")

    if solver.StatusName(status) == "INFEASIBLE":
        print(f"Model is infeasible; {solver.SufficientAssumptionsForInfeasibility()}")

# %%
