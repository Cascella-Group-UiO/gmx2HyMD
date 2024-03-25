import numpy as np


class GroAtom:
    def __init__(
        self, resid, resname, atom_name, index, x, y, z, vx=0.0, vy=0.0, vz=0.0
    ):
        self.resid = resid
        self.resname = resname
        self.atom_name = atom_name
        self.index = index
        self.x = x
        self.y = y
        self.z = z
        self.vx = vx
        self.vy = vy
        self.vz = vz

    @classmethod
    def parse_line(cls, line):
        line_length = len(line)
        resid = int(line[:5])
        resname = line[5:10].strip()
        atom_name = line[10:15].strip()
        index = int(line[15:20])
        x = float(line[20:28])
        y = float(line[28:36])
        z = float(line[36:44])
        vx = float(line[44:52]) if line_length == 69 else 0.0
        vy = float(line[52:60]) if line_length == 69 else 0.0
        vz = float(line[60:68]) if line_length == 69 else 0.0
        if line_length not in (45, 69):
            raise ValueError(
                f"Gro file line not formatted correctly:\n"
                f"{line}"
                "The line lenght is {line_length} "
                "while it should be 45 for a gro file containing positions only "
                "or 69 for a gro file containing both positions and velocities."
            )
        return cls(resid, resname, atom_name, index, x, y, z, vx, vy, vz)


def load_gro(filename: str) -> tuple[list[GroAtom], np.ndarray]:
    """Parse gro file"""
    with open(filename, "r") as infile:
        lines = infile.readlines()

    atom_list = []
    box_size = np.array(lines[-1].split(), dtype=float)
    for line in lines[2:-1]:
        atom_list.append(GroAtom.parse_line(line))
    print(f"GRO file {filename} loaded... ")
    return atom_list, box_size
