import h5py
import numpy as np

from .gro_utils import GroAtom
from .pdb_utils import PdbAtom
from .system_properties import System

AtomList = list[GroAtom] | list[PdbAtom]
# Protein beads
one_to_three = {
    # fmt: off
    "A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE",
    "G": "GLY", "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU",
    "M": "MET", "N": "ASN", "P": "PRO", "Q": "GLN", "R": "ARG",
    "S": "SER", "T": "THR", "V": "VAL", "W": "TRP", "Y": "TYR",
}
three_to_one = {one_to_three[k]: k for k in one_to_three}


def write_coordinates(
    atoms: AtomList,
    system: System,
    box: np.ndarray,
    electric_label: bool,
    outfile: str,
):
    with h5py.File(outfile, "w") as f_h5:
        f_h5.attrs["box"] = box[:3]  # Only cubic boxes for now
        f_h5.attrs["n_molecules"] = system.n_mol

        n_atoms = len(atoms)
        f_h5.create_dataset("indices", data=np.arange(n_atoms), dtype="i")
        f_h5.create_dataset("names", data=system.names, dtype="S5")
        f_h5.create_dataset("types", data=system.types, dtype="i")
        f_h5.create_dataset("molecules", data=system.molecules, dtype="i")
        f_h5.create_dataset("masses", data=system.masses, dtype="float32")
        write_pos_and_vel(f_h5, n_atoms, atoms)

        if electric_label:
            f_h5.create_dataset("charge", data=system.charges, dtype="float32")

        # TODO: add bonds dataset back? Useless if using index based topology
        # MAX_N_BONDS = 4  # can potentially be changed to a higher number if needed
        # dset_bonds = f_h5.create_dataset("bonds", (n_atoms, MAX_N_BONDS), dtype="i")
        # dset_bonds[...] = access_top_molecule_itps_gen_whole_atomBondIdx(
        #     molecule_list, atomcsv
        # )


def write_pos_and_vel(f_h5: h5py.File, n_atoms: int, atoms: AtomList):
    dset_pos = f_h5.create_dataset("coordinates", (1, n_atoms, 3), dtype="float64")
    dset_vel = f_h5.create_dataset("velocities", (1, n_atoms, 3), dtype="float64")
    for idx, atom in enumerate(atoms):
        dset_pos[0, idx, :] = np.array([atom.x, atom.y, atom.z])
        dset_vel[0, idx, :] = np.array([atom.vx, atom.vy, atom.vz])
