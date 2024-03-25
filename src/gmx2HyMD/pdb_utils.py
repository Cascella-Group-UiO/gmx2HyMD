from typing import Optional

import numpy as np


class PdbAtom:
    def __init__(
        # fmt: off
        # these argmunets are all 'str'
        self, label, index, name, residue, resid,
        x, y, z, occu, temp, element,
    ):
        self.label = str(label)
        self.index = int(index)
        self.name = str(name)
        self.residue = str(residue)
        self.resid = int(resid)
        self.x = float(x) / 10.0
        self.y = float(y) / 10.0
        self.z = float(z) / 10.0
        self.vx = 0.0
        self.vy = 0.0
        self.vz = 0.0
        self.occu = float(occu)
        self.temp = float(temp)
        self.element = str(element)


def load_pdb(
    filename: str, cli_box: Optional[np.ndarray]
) -> tuple[list[PdbAtom], np.ndarray]:
    atom_list = []
    box = None
    with open(filename, "r") as infile:
        for line in infile:
            if line.startswith("END"):
                break
            if line.startswith("REMARK"):
                continue
            if line.startswith("CRYST1"):
                box = np.array(line.split()[1:4], dtype=float) / 10.0
                continue
            if line.startswith("TER"):
                continue
            args = line.split()
            if len(args) != 11:
                raise ValueError(
                    "Only pdb files with exactly 11 entries per atom line are allowed:\n"
                    "label, index, atom name, residue name, residue index, "
                    "pos_x, pos_y, pos_z, occupation, beta factor, element"
                )
            atom_list.append(PdbAtom(*args))

    if box is None:
        if cli_box is not None:
            box = np.array(cli_box)
        else:
            raise ValueError(
                "Trying to convert a pdb file, but box info was either not available in the file "
                "('CRYST1' line missing) or was not provided via command line."
            )
    print(f"PDB file {filename} loaded... ")
    return atom_list, box
