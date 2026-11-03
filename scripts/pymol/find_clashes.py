"""Find exceptionally close heavy-atom pairs in a loaded PyMOL structure.

Load this file in PyMOL with ``run path/to/find_clashes.py``. It registers the
``find_clashes`` command. The search uses PyMOL's native ``cmd.find_pairs``
routine, so it checks nearby atom pairs instead of comparing every atom with
every other atom in Python.
"""

from __future__ import print_function

import csv
import math

from pymol import cmd


CSV_FIELDS = (
    "model_1", "atom_index_1", "segment_1", "chain_1", "residue_id_1",
    "residue_name_1", "atom_name_1", "element_1",
    "model_2", "atom_index_2", "segment_2", "chain_2", "residue_id_2",
    "residue_name_2", "atom_name_2", "element_2", "distance_angstrom",
)


def _as_bool(value):
    """Accept Python booleans, integers, and common command-line spellings."""
    return str(value).strip().lower() not in ("0", "false", "no", "off", "")


def _distance(atom_a, atom_b):
    """Return the Euclidean distance between two atom records."""
    dx = atom_a[7] - atom_b[7]
    dy = atom_a[8] - atom_b[8]
    dz = atom_a[9] - atom_b[9]
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def find_clashes(selection="all", cutoff=1.0, state=1,
                 object_name="clash_lines", exclude_same_residue=1,
                 include_hydrogens=0, csv_file=""):
    """Draw and report atom pairs closer than ``cutoff`` angstroms.

    Args:
        selection: PyMOL selection to scan. Defaults to all loaded atoms.
        cutoff: Maximum atom-center distance to report, in angstroms.
        state: PyMOL state to scan. Defaults to state 1.
        object_name: Name of the red distance object created in PyMOL.
        exclude_same_residue: Set to 1 to omit pairs from the same residue.
        include_hydrogens: Set to 1 to include hydrogen atoms.
        csv_file: Optional path for a CSV report. Empty means no file is written.

    Returns:
        Number of close atom pairs found.
    """
    try:
        cutoff = float(cutoff)
        state = int(state)
    except (TypeError, ValueError):
        raise ValueError("cutoff must be a number and state must be an integer")

    selection = str(selection).strip()
    object_name = str(object_name).strip()
    if not selection:
        raise ValueError("selection cannot be empty")
    if cutoff <= 0:
        raise ValueError("cutoff must be greater than zero")
    if state < 1:
        raise ValueError("state must be 1 or greater")
    if not object_name:
        raise ValueError("object_name cannot be empty")

    exclude_same_residue = _as_bool(exclude_same_residue)
    include_hydrogens = _as_bool(include_hydrogens)

    atom_selection = "({})".format(selection)
    if not include_hydrogens:
        atom_selection += " and not elem H"

    # Keep metadata and coordinates in one lookup, keyed by PyMOL model/index.
    # get_coords returns coordinates in the same order as iterate and accounts
    # for object transforms, which iterate_state coordinates do not.
    atom_refs = []
    atom_info = {}
    cmd.iterate(
        atom_selection,
        "atom_refs.append((str(model), int(index))); "
        "atom_info[(str(model), int(index))] = "
        "(segi, chain, resi, resn, name, elem, alt)",
        space={"atom_refs": atom_refs, "atom_info": atom_info},
    )

    atom_data = {}
    if atom_refs:
        coordinates = cmd.get_coords(atom_selection, state)
        if coordinates is None or len(coordinates) != len(atom_refs):
            raise RuntimeError(
                "PyMOL returned a different number of coordinates and atoms; "
                "try scanning one object at a time."
            )
        for atom_ref, coordinate in zip(atom_refs, coordinates):
            atom_data[atom_ref] = atom_info[atom_ref] + (
                float(coordinate[0]),
                float(coordinate[1]),
                float(coordinate[2]),
            )

    candidate_pairs = []
    if atom_data:
        candidate_pairs = cmd.find_pairs(
            atom_selection,
            atom_selection,
            state1=state,
            state2=state,
            cutoff=cutoff,
            mode=0,
        )

    clashes = []
    seen_pairs = set()

    for pair in candidate_pairs:
        ref_a = (str(pair[0][0]), int(pair[0][1]))
        ref_b = (str(pair[1][0]), int(pair[1][1]))
        if ref_a == ref_b:
            continue

        pair_key = tuple(sorted((ref_a, ref_b)))
        if pair_key in seen_pairs:
            continue
        seen_pairs.add(pair_key)

        if ref_a not in atom_data or ref_b not in atom_data:
            continue

        atom_a = atom_data[ref_a]
        atom_b = atom_data[ref_b]

        # The residue key includes model, segment, chain, residue ID, and name.
        residue_a = (ref_a[0],) + atom_a[:4]
        residue_b = (ref_b[0],) + atom_b[:4]
        if exclude_same_residue and residue_a == residue_b:
            continue

        distance = _distance(atom_a, atom_b)
        if distance >= cutoff:
            continue

        clash = {
            "model_1": ref_a[0],
            "atom_index_1": ref_a[1],
            "segment_1": atom_a[0],
            "chain_1": atom_a[1],
            "residue_id_1": atom_a[2],
            "residue_name_1": atom_a[3],
            "atom_name_1": atom_a[4],
            "element_1": atom_a[5],
            "model_2": ref_b[0],
            "atom_index_2": ref_b[1],
            "segment_2": atom_b[0],
            "chain_2": atom_b[1],
            "residue_id_2": atom_b[2],
            "residue_name_2": atom_b[3],
            "atom_name_2": atom_b[4],
            "element_2": atom_b[5],
            "distance_angstrom": round(distance, 3),
        }
        clashes.append(clash)

    clashes.sort(key=lambda item: item["distance_angstrom"])

    atom_selection_name = object_name + "_atoms"
    residue_selection_name = object_name + "_residues"

    # Keep the view focused on this run: show sticks only for clash residues.
    cmd.hide("labels", "all")
    cmd.hide("sticks", "all")
    cmd.delete(object_name + "_context")
    cmd.delete(atom_selection_name)
    cmd.delete(residue_selection_name)

    # Replace the output object from a previous run, if present.
    cmd.delete(object_name)

    if clashes:
        atom_selections = set()
        for clash in clashes:
            selection_1 = "model {} and index {}".format(
                clash["model_1"], clash["atom_index_1"]
            )
            selection_2 = "model {} and index {}".format(
                clash["model_2"], clash["atom_index_2"]
            )
            atom_selections.add("({})".format(selection_1))
            atom_selections.add("({})".format(selection_2))
            cmd.distance(object_name, selection_1, selection_2)

        cmd.select(atom_selection_name, " or ".join(sorted(atom_selections)))
        cmd.select(
            residue_selection_name,
            "byres ({})".format(atom_selection_name),
        )
        sidechain_selection = "({}) and sidechain and not elem H".format(
            residue_selection_name
        )
        # Also show a clashing atom if it is part of the backbone.
        stick_selection = "({}) or ({})".format(
            sidechain_selection, atom_selection_name
        )
        cmd.show("sticks", stick_selection)
        cmd.color(
            "orange",
            "({}) and not elem H".format(residue_selection_name),
        )
        cmd.color("yellow", atom_selection_name)

        # Match the original, easy-to-see PyMOL measurement style.
        cmd.hide("labels", object_name)
        cmd.color("red", object_name)
        cmd.set("dash_width", 4)
        cmd.set("dash_gap", 0)
        cmd.zoom(object_name)

    if str(csv_file).strip():
        with open(str(csv_file).strip(), "w", newline="") as output_handle:
            writer = csv.DictWriter(output_handle, fieldnames=CSV_FIELDS)
            writer.writeheader()
            for clash in clashes:
                writer.writerow({field: clash[field] for field in CSV_FIELDS})

    print(
        "Detected {} atom pair(s) closer than {:.3f} angstroms "
        "in state {} (selection: {}).".format(len(clashes), cutoff, state, selection)
    )
    if not clashes:
        print("No clashes found at the chosen cutoff.")
    if str(csv_file).strip():
        print("Wrote clash details to: {}".format(str(csv_file).strip()))

    return len(clashes)


cmd.extend("find_clashes", find_clashes)
