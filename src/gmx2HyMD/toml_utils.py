import os

import tomlkit


def toml_add_lines(input_list):
    """"""
    out_arr = tomlkit.array()
    for i in input_list:
        out_arr.add_line(i)
    out_arr.add_line(indent="")
    return out_arr


def add_mol_section(section, params, name):
    if params:
        converted_list = [list(elem.__dict__.values()) for elem in params]
        section.add(name, toml_add_lines(converted_list))
        section.add(tomlkit.nl())


def write_itp_toml(molname, params, dirname):
    doc = tomlkit.document()
    molecule = tomlkit.table()
    molecule.add("atomnum", params.n_atoms)

    # TODO: Ideally we would iterate over the class
    add_mol_section(molecule, params.atoms, "atoms")
    add_mol_section(molecule, params.bonds, "bonds")
    add_mol_section(molecule, params.angles, "angles")
    add_mol_section(molecule, params.dihedrals, "dihedrals")

    doc.add(molname, molecule)
    with open(f"{dirname}/{molname}.toml", "w") as outfile:
        tomlkit.dump(doc, outfile)


def write_topol_toml(molecule_list, out_name):
    doc = tomlkit.document()
    system = tomlkit.table()
    system.add("molecules", toml_add_lines([[k, v] for k, v in molecule_list.items()]))
    system.add(
        "include", toml_add_lines([f"./{key}.toml" for key in molecule_list.keys()])
    )
    doc.add("system", system)

    with open(out_name, "w") as outfile:
        tomlkit.dump(doc, outfile)


def write_topology(molecule_list, topol, outfile):
    dirname = os.path.dirname(os.path.abspath(outfile))
    for molecule in molecule_list:
        try:
            # write single file toml for each molecule type
            write_itp_toml(molecule, topol[molecule], dirname)
        except ValueError:
            raise ValueError(
                f"{molecule}: missing itp file where the {molecule} moleculetype is defined!"
            )
    # write topol.toml
    write_topol_toml(molecule_list, outfile)
