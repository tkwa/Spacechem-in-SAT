

### Features not implemented yet

- Input and output.
    - Write define_input_pattern()
    - Enforce atoms only change upon activation
- Bond+ and Bond- commands.
- Rotation.
- Multiple bonds and bond limits.
- Some way to import puzzles and export solutions.
- Two waldos. (It would be simpler to add the extra waldo if it could never grab or output.)
- Flip flops.
- Sensors.
- Sync
- Fuse and split.
- Teleporter
- Random inputs

# Pending refactor

It would potentially be good to move from object-oriented (waldo, atom, cell) to grouping commands by functionality.
 - waldo pathing
 - bonding (restrictions on which bonds form when)
 - molecule movement.

### Limitations (not mentioned in "Features not implemented yet")

- Input and output cannot happen on the same cycle