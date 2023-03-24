
A solver for SpaceChem research puzzles using ORTools.

## How to use

Don't, it's not functional yet. If you insist, look at the tests in `test.py` and write your own puzzle spec.

## Features

- Waldos can move around, grab and drop atoms.
- Atoms can be connected by single bonds.
- Input and output.

## Features not implemented yet

- Bond+ and Bond- commands.
- Some way to import puzzles and export solutions.
- Loop checking (for finding the lowest cycle solution).
- Input beta and output omega zones.
- Rotation.
- Multiple bonds and bond limits.
- Output isomorphisms.
- Collisions during rotation.
- Two waldos. (It would be simpler to add the extra waldo if it could never grab.)
- Flip flops.
- Sensors.
- Sync
- Fuse and split.
- Teleporter
- Random inputs

## Pending refactor

It would potentially be good to move from object-oriented (waldo, atom, cell) to grouping commands by functionality.
 - waldo pathing
 - bonding (restrictions on which bonds form when)
 - molecule movement.

## Limitations (not mentioned in "Features not implemented yet")

- Input and output cannot happen on the same cycle. This is because input is implemented by allocating new atoms, and output is implemented by deallocating atoms.
- Output is lenient; when there are multiple molecules in the output zone, SpaceChem sometimes outputs the wrong molecule