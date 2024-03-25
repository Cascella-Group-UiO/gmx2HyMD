from dataclasses import dataclass

import numpy as np

from .topol_utils import Topol


@dataclass
class System:
    n_mol: int
    names: np.ndarray
    types: np.ndarray
    molecules: np.ndarray
    masses: np.ndarray
    charges: np.ndarray


def get_system_properties(
    molecule_list: dict[str, int], topol: dict[str, Topol]
) -> System:
    # TODO: read types from database given names otherwise create name to type converter
    # TODO: introduce flag to use given forcefield (useful for proteins, lipids, etc.)
    names, molecules, masses, charges = [], [], [], []
    exclude = ("W", "SOL", "ION", "HOH")
    n_mol = 0
    for molname, molnum in molecule_list.items():
        # type(topol[molname]) == Topol, type(atom) == list(ItpAtom)
        names += molnum * [
            atom.atomtype if atom.resname not in exclude else atom.atomname
            for atom in topol[molname].atoms
        ]

        for i in range(molnum):
            molecules += [(n_mol + i)] * topol[molname].n_atoms
        atom_masses = [atom.mass for atom in topol[molname].atoms]
        masses += atom_masses * molnum
        n_mol += molnum

        charges += molnum * [atom.charge for atom in topol[molname].atoms]

    names = np.array(names, dtype="S5")
    molecules = np.array(molecules)
    masses = np.array(masses)
    charges = np.array(charges)

    # types
    _, idx = np.unique(names, return_index=True)
    unique_names = names[np.sort(idx)]
    name_to_type = {name: t for t, name in enumerate(unique_names)}
    types = np.array([name_to_type[name] for name in names])

    return System(n_mol, names, types, molecules, masses, charges)
