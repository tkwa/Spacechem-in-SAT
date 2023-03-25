
A solver for SpaceChem research puzzles using ORTools.

## How to use

Importing puzzles is not supported yet. Run `it-takes-three.py` to generate a solution to an easy puzzle; it should take about 100 seconds.

## Features

- Waldos can move around, grab and drop atoms.
- Atoms can be connected by single bonds.
- Input and output.
- Bond+ and Bond- commands.

## Features not implemented yet

- Some way to import puzzles and export solutions.
- Output omega zone.
- Rotation.
- Multiple bonds and bond limits.
- Output isomorphisms. Currently not even translation of an output molecule is supported.
- Collisions during rotation.
- Two waldos. (It would be simpler to add the extra waldo if it could never grab.)
- Flip flops.
- Sensors.
- Sync.
- Fuse and split.
- Teleporter.
- Random inputs.

## Pending refactor

It would potentially be good to move from object-oriented (waldo, atom, cell) to grouping commands by functionality.
 - waldo pathing
 - bonding (restrictions on which bonds form when)
 - molecule movement.

## Limitations (not mentioned in "Features not implemented yet")

- Input and output cannot happen on the same cycle. This is because input is implemented by allocating new atoms, and output is implemented by deallocating atoms.
- Output is lenient; when there are multiple molecules in the output zone, SpaceChem sometimes outputs the wrong molecule. The solver merely requires that the correct molecule is in the output zone.
- The solver cannot easily translate solutions across time (and presumably space).