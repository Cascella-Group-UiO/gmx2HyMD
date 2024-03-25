import os
from typing import Callable, TypeVar

import parse


class Topol:
    def __init__(self, atoms, bonds=None, angles=None, dihedrals=None):
        if bonds is None:
            bonds = []
        if angles is None:
            angles = []
        if dihedrals is None:
            dihedrals = []

        self.n_atoms = len(atoms)
        self.atoms = atoms
        self.bonds = bonds
        self.angles = angles
        self.dihedrals = dihedrals


class ItpAtom:
    def __init__(
        self, index, atomtype, resnr, resname, atomname, cgnr, charge, mass=None
    ):
        self.index = int(index)
        self.atomtype = atomtype
        self.resnr = int(resnr)
        self.resname = resname
        self.atomname = atomname
        self.cgnr = int(cgnr)
        self.charge = float(charge)

        if mass is not None:
            self.mass = float(mass)
            return

        # MARTINI 2
        self.mass = 72.0

        # MARTINI 3
        # if atomtype.startswith("T"):
        #     self.mass = 32.0
        # elif atomtype.startswith("S"):
        #     self.mass = 54.0
        # else:
        #     self.mass = 72.0


class ItpBond:
    def __init__(self, atom1, atom2, func, lenght=999, strength=999):
        self.atom1 = int(atom1)
        self.atom2 = int(atom2)
        self.func = int(func)
        self.length = float(lenght)
        self.strength = float(strength)


class ItpAngle:
    def __init__(self, atom1, atom2, atom3, func, lenght=999, strength=999):
        self.atom1 = int(atom1)
        self.atom2 = int(atom2)
        self.atom3 = int(atom3)
        self.func = int(func)
        self.length = float(lenght)
        self.strength = float(strength)


class ItpDihedral:
    def __init__(self, atom1, atom2, atom3, atom4, func, lenght=999, strength=999):
        self.atom1 = int(atom1)
        self.atom2 = int(atom2)
        self.atom3 = int(atom3)
        self.atom4 = int(atom4)
        self.func = int(func)

        # Improper dihedral
        if self.func == 2:
            self.length = float(lenght)
            self.strength = float(strength)


ItpSection = TypeVar("ItpSection", int, ItpAtom, ItpBond, ItpAngle, ItpDihedral)


def load_itp_section(data: list[str], class_or_fun: Callable) -> list[ItpSection]:
    """ """
    output_list = []
    ifdef = False
    inner_ifdef = False
    for line in data:
        # CHECK: what's the best way to order these if-statements?
        # break if next section starts
        if line.startswith("["):
            break

        # Skip #ifdef blocks
        if line.startswith("#if"):
            if ifdef:
                inner_ifdef = True
            else:
                ifdef = True
            continue
        if line.startswith("#endif"):
            if inner_ifdef:
                continue
            ifdef = False
            continue
        if ifdef:
            continue

        # Skip comments or empty lines
        if line in ("\n", " \n") or line.startswith(";"):
            continue
        demoline = line.split(";")[0].split()
        if not demoline:
            continue

        # CHECK: I feel like classes only add complexity without benifits here
        output_list.append(class_or_fun(*demoline))
    return output_list


def parse_itp(
    itp_file: str, molecule_list: dict[str, int], elec_label: bool
) -> tuple[dict[str, tuple[ItpAtom, ItpBond, ItpAngle, ItpDihedral]], bool]:
    # tuple[dict[str, tuple[]]] :
    atoms_list, bonds_list, angles_list, dih_list = [], [], [], []

    print("Reading ITP file:", itp_file)
    with open(itp_file, "r") as infile:
        lines = infile.readlines()

    moltype_idx = []
    for i, line in enumerate(lines):
        if line.startswith("[ moleculetype") or line.startswith("[moleculetype"):
            moltype_idx.append(i + 1)

    if not moltype_idx:
        raise ValueError("Missing [ moleculetype ] section in {itp_file}.")

    molecules = {}
    for i, idx in enumerate(moltype_idx):
        try:
            end = moltype_idx[i + 1]
        except IndexError:
            end = None
        sel_lines = lines[idx:end]

        molname = load_itp_section(sel_lines, lambda x, _: x)[0]
        if molname in molecule_list:
            sections = {}
            for j, line in enumerate(sel_lines):
                pattern = parse.parse("[{}]\n", line)
                if pattern is not None:
                    sections[pattern[0].strip()] = j + 1

            if "atoms" in sections:
                start = sections["atoms"]
                atoms_list = load_itp_section(sel_lines[start:], ItpAtom)
                if any([atom.charge for atom in atoms_list]):
                    elec_label = True
            if "bonds" in sections:
                start = sections["bonds"]
                bonds_list = load_itp_section(sel_lines[start:], ItpBond)
            if "angles" in sections:
                start = sections["angles"]
                angles_list = load_itp_section(sel_lines[start:], ItpAngle)
            if "dihedrals" in sections:
                start = sections["dihedrals"]
                dih_list = load_itp_section(sel_lines[start:], ItpDihedral)
            molecules[molname] = (atoms_list, bonds_list, angles_list, dih_list)
    return molecules, elec_label


def load_top_params(
    molecule_list: dict[str, int], itp_paths: list[str], elec_label: bool
) -> tuple[dict[str, Topol], int, bool]:
    topol = {}
    for itp in itp_paths:
        # FIXME: if file with LJ params do something else?
        # if os.path.basename(itp) in [
        # "martini.itp",
        # "martini_v2.2.itp",
        params, elec_label = parse_itp(itp, molecule_list, elec_label)
        for molname in params:
            topol[molname] = Topol(*params[molname])

    topol_atoms = 0
    print("System composition:")
    for molname in molecule_list:
        n_mol = molecule_list[molname]
        # if molname in ["W", "NA", "CL"] and molname not in topol:
        #     topol[molname] = Topol(*single_bead_itp(molname))
        topol_atoms += n_mol * topol[molname].n_atoms
        print(f"    {molname:10}\t{n_mol:>6}")
    return topol, topol_atoms, elec_label


def load_top(filename: str) -> tuple[dict[str, int], list[str]]:
    """
    Lookup top file's [ molecules ] section and included itps
    """
    itp_molecules = {}
    itp_paths = []
    with open(filename, "r") as infile:
        lines = infile.readlines()

    try:
        idx = lines.index("[ molecules ]\n")
    except ValueError:
        try:
            idx = lines.index("[molecules]\n")
        except:
            raise ValueError(f"[ molecules ] section missing in {filename}.")

    for line in lines[:idx]:
        if line.startswith("#"):
            line_wo_comment = line.split(";")[0]
            path = parse.parse("#include {}\n", line_wo_comment)
            if path != None:
                itp_paths.append(
                    f"{os.path.dirname(os.path.abspath(filename))}/{path[0][1:-1]}"
                )  # strip sorrounding quotes

    for line in lines[idx + 1 :]:
        if line == "\n":
            # cut off at the first empty line
            break
        elif line.startswith(";"):
            # Skip comments
            continue
        else:
            # make sure only the first two elements are used
            molname, molnum = line.split()[:2]
            if molname in itp_molecules:
                itp_molecules[molname] += int(molnum)
            else:
                itp_molecules[molname] = int(molnum)
    print(f"TOP file {filename} loaded... ")
    return itp_molecules, itp_paths
