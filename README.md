# PyMOL Structure Clash Checker

A small PyMOL script that finds exceptionally close atom pairs in a loaded molecular structure, colors directly involved residues orange and clashing atoms yellow, shows their side chains as sticks, and draws red PyMOL distance measurements. It can also save atom identities and distances to CSV.

## What it checks

By default, the command scans all atoms in the PyMOL session, excludes hydrogen atoms, skips pairs within the same residue, and reports pairs whose atom centers are less than **1.0 angstrom apart**. PyMOL's native neighbor-search routine supplies nearby candidates, so the script does not compare every atom pair in Python.

This is a screen for exceptionally short contacts. It does not calculate van der Waals overlap, validate bond geometry, or replace a complete structure-validation workflow. The cutoff is an absolute distance in angstroms. If you raise it substantially, covalently bonded atoms in neighboring residues may also be reported; the script only excludes same-residue pairs, not explicit bonds.

## Requirements

- PyMOL with its Python API available
- No additional Python packages
- A structure already loaded in PyMOL (for example, a PDB or mmCIF file)

## Quick start

1. Load a structure in PyMOL, for example:

   ```text
   load path/to/structure.pdb, my_structure
   ```

   You can also load an mmCIF file or use any structure already open in the session.

2. From the repository root, load the script in the PyMOL command line:

   ```text
   run scripts/pymol/find_clashes.py
   ```

3. Run the default scan:

   ```text
   find_clashes
   ```

   PyMOL prints the number of pairs found and draws them in a red distance-measurement object named `clash_lines`. The directly involved residues are colored orange, their clashing atoms are colored yellow, and only their side chains are shown as sticks. Labels are hidden, other stick representations are hidden, and the view zooms to the clash lines.

To run the same quick-start commands from a PyMOL script after loading a structure, use:

```text
@examples/check_loaded_structure.pml
```

## Command options

The registered command accepts these positional arguments in this order:

```text
find_clashes [selection [, cutoff [, state [, object_name [, exclude_same_residue [, include_hydrogens [, csv_file]]]]]]]
```

| Argument | Default | Meaning |
| --- | --- | --- |
| `selection` | `all` | Any valid PyMOL atom selection. For example, `polymer.protein`, `chain A`, or a named object. |
| `cutoff` | `1.0` | Maximum atom-center distance to report, in angstroms. Only distances strictly below the cutoff are included. |
| `state` | `1` | PyMOL state to scan. |
| `object_name` | `clash_lines` | Name of the red PyMOL distance object. A previous object with this name is replaced. |
| `exclude_same_residue` | `1` | Use `1` to skip atom pairs in the same model, segment, chain, residue ID, and residue name; use `0` to include them. |
| `include_hydrogens` | `0` | Use `1` to include hydrogen atoms; by default, atoms with element H are excluded. |
| `csv_file` | empty | Optional output path for a CSV report. Relative paths are resolved from PyMOL's current working directory. |

Examples:

```text
# Scan just the protein atoms at the default cutoff.
find_clashes polymer.protein, 1.0

# Scan chain A, save the lines under a different name, and write a CSV report.
find_clashes (chain A), 1.0, 1, chain_A_clashes, 1, 0, chain_A_clashes.csv

# Include same-residue pairs and hydrogen atoms.
find_clashes all, 1.0, 1, all_atom_clashes, 0, 1

```

The CSV contains both atoms' model, atom index, segment, chain, residue, atom name, and element, plus their separation in angstroms. It is useful when you need to inspect or sort the reported pairs outside PyMOL.

## Notes and assumptions

- The default selection is `all`, so non-protein atoms in the session (such as ligands, waters, or ions) are included. Use `polymer.protein` or another PyMOL selection to narrow the clash scan.
- Atoms in alternate conformations are included if they match the selection. If the structure has overlapping alternate locations, choose one conformer in the selection before running the scan.
- For a multi-state object, only the requested state is checked; the default is state 1.
- Same-residue exclusion is based on model, segment, chain, residue ID, and residue name. It does not infer chemical connectivity.
- The script colors heavy atoms in each directly involved residue orange and the clashing atoms yellow. It shows the residues' side chains as sticks; a clashing backbone atom is also shown as a stick so the measured pair is visible. The clash atom and residue selections are named `<object_name>_atoms` and `<object_name>_residues` (by default, `clash_lines_atoms` and `clash_lines_residues`).
- To keep the view focused, the script hides stick representations and labels across the session before drawing the current clash results. Other stick representations are therefore hidden when the command runs.
- The output object name is replaced each time the command runs. Choose another `object_name` if you want to keep separate scan results in the same session.

## Project files

```text
pymol-structure-clash-checker/
├── README.md
├── .gitignore
├── examples/
│   └── check_loaded_structure.pml
└── scripts/
    └── pymol/
        └── find_clashes.py
```

No structure files or generated results are included, so you can add the project to a public repository without including your own input structures.

## PyMOL references

- [PyMOL `find_pairs` API](https://wiki.pymol.org/index.php/Find_pairs)
- [PyMOL `get_coords` API](https://wiki.pymol.org/index.php/Get_coords)
- [PyMOL command reference](https://pymol.org/pymol-command-ref.html)
