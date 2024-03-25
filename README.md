# gmx2HyMD
gmx2HyMD is a CLI utility that converts
[MARTINI](http://www.cgmartini.nl/index.php/martini) coordinates and topologies
(obtained for example with [CHARMM-GUI](https://www.charmm-gui.org/)) to inputs
for [HyMD](https://github.com/Cascella-Group-UiO/HyMD) and [∂-HyMD](https://github.com/Cascella-Group-UiO/Diff-HyMD).

After cloning the repository, the program can be installed with pip.\
Optionally you can also install it in a virtual environment
```terminal
python -m venv <your_venv_dir> --upgrade-deps
source <your_venv_dir>/bin/activate
```
Then
```terminal
cd gmx2HyMD
pip install .
```

You should now have access to the `gmx2HyMD` command
```terminal
gmx2HyMD -f <input>.gro -p <topol>.top
```

Check the the other available flags with
```terminal
gmx2HyMD --help
```
In the `test_system` directory you can find example input files to use with `gmx2HyMD`.
