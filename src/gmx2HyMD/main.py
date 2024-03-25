import argparse
import os
import time

import numpy as np

from .gro_utils import load_gro
from .h5_utils import three_to_one, write_coordinates
from .pdb_utils import load_pdb
from .system_properties import get_system_properties
from .toml_utils import write_topology
from .topol_utils import load_top, load_top_params


def write_simulation_parameters(names: np.ndarray):
    template = os.path.abspath(os.path.dirname(__file__)) + "/template.toml"
    out_lines = []

    _, idx = np.unique(names, return_index=True)
    unique_names = np.array(names[np.sort(idx)], dtype=str)

    i, j = np.triu_indices(unique_names.size, k=1)
    chi_pairs = np.column_stack(
        (
            [f"'{n}'" for n in unique_names[i]],
            [f"'{n}'" for n in unique_names[j]],
        )
    )
    filter = np.isin(unique_names, ("W", "NA", "NA+", "CL", "CL-", "CA"))

    with open(template, "r") as infile:
        for line in infile:
            if line.startswith("thermostat_coupling_groups"):
                new_line = (
                    "thermostat_coupling_groups = [\n"
                    f"    {np.array2string(unique_names[~filter], separator=', ')},\n"
                    f"    {np.array2string(unique_names[filter], separator=', ')},\n]\n"
                )
                out_lines.append(new_line)
            # chi matrix
            elif line.startswith("]"):
                # TODO: replace 10 with values from database
                new_line = "".join(
                    [f"    [{i:>7}, {j:>7}, {10:>5}],\n" for i, j in chi_pairs]
                )
                out_lines.append(new_line + "]\n")
            else:
                out_lines.append(line)

    with open("options.toml", "w") as outfile:
        outfile.write("".join(out_lines))


def user_input() -> argparse.Namespace:
    description = "Convert GROMACS coordinates (gro/pdb) and topology (top) to HyMD compatible H5 and toml files."
    ap = argparse.ArgumentParser(description=description)
    ap.add_argument(
        "-f",
        "--file",
        required=True,
        dest="input",
        type=str,
        help="Input .gro/.pdb file",
    )
    ap.add_argument(
        "-p", "--top", required=True, dest="top", type=str, help="Input .top file"
    )
    ap.add_argument(
        "-b",
        "--box",
        nargs=3,
        type=float,
        default=None,
        help="Box size, takes 3 inputs (x, y, z). Required if the input pdb file does not provide it.",
    )
    ap.add_argument(
        "-oc",
        dest="out_h5",
        default=None,
        type=str,
        help="Output H5MD file (defaults to 'output.h5')",
    )
    ap.add_argument(
        "-op",
        dest="out_toml",
        default=None,
        type=str,
        help="Output toml file (defaults to 'topol.toml')",
    )
    ap.add_argument(
        "--nopar",
        action="store_true",
        dest="no_params",
        help="Disable 'options.toml' file output.",
    )
    ap.add_argument(
        "--notop",
        action="store_true",
        dest="no_topol",
        help="Disable 'topol.toml' file output.",
    )
    return ap.parse_args()


def main():
    """Read coarse grained .gro/pdb and .top files from gmx Martini and write HyMD input (.h5) and topology (.toml)."""

    start_time = time.time()
    args = user_input()
    elec_label = False  # TODO: use zero charges if no charges are provided instead and parse correctly in HyMD

    # Get topology parameters
    molecule_list, itp_paths = load_top(args.top)
    topol, topol_atoms, elec_label = load_top_params(
        molecule_list, itp_paths, elec_label
    )
    molecule_idx = {mol: idx for idx, mol in enumerate(molecule_list)}

    # Get atom positions
    base, ext = os.path.splitext(args.input)
    if ext == ".gro":
        atoms, box = load_gro(args.input)
        protein_in_top = any([atom.resname in three_to_one for atom in atoms])
        if protein_in_top:
            print(
                "WARNING: converting system containg protein molecules. Atoms will have the same ordering as in the input coordinate file. "
                "Depending on the system, this might lead to inconsistencies with the output topology."
            )
        else:
            atoms.sort(key=lambda x: molecule_idx[x.resname])
    elif ext == ".pdb":
        atoms, box = load_pdb(args.input, args.box)
        protein_in_top = any([atom.residue in three_to_one for atom in atoms])
        if protein_in_top:
            print(
                "WARNING: converting system containg protein molecules. Atoms will have the same ordering as in the input coordinate file. "
                "Depending on the system, this might lead to inconsistencies with the output topology."
            )
        else:
            atoms.sort(key=lambda x: molecule_idx[x.residue])
    else:
        raise ValueError(
            f"Input coordinate file extension should either be .gro or .pdb. Got {ext}."
        )

    # Check that the number of atoms is the same in both files
    if topol_atoms != len(atoms):
        raise ValueError(
            f"Atom number mismatch in {args.input} (# {len(atoms)}) and {args.top} (# {topol_atoms})."
        )

    # write output files
    if args.out_h5 is None:
        args.out_h5 = "./input.h5"
    if args.out_toml is None:
        args.out_toml = f"./topol.toml"

    system = get_system_properties(molecule_list, topol)
    # properties = (n_mol, names, types, molecules, masses, charges)

    if not args.no_topol:
        write_topology(molecule_list, topol, args.out_toml)
    if not args.no_params:
        write_simulation_parameters(system.names)

    write_coordinates(atoms, system, box, elec_label, args.out_h5)
    print(f"Conversion took {time.time() - start_time} s.")


if __name__ == "__main__":
    main()
