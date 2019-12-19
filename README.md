# CIF2ZTB: Tabulated potential of material structures

A CUDA program to calculate tabulated potentials for adsorption simulation using an implicit rigid adsorbent material defined by its potential energy surface. Based on the Fortran code of [MCCCS-MN](https://github.com/SiepmannGroup/VLE-Validation/tree/master/MCCCS-MN) and refactored in C++/CUDA.

## Compilation and usage
Use `make` to compile the program. 

Usage: `cif2ztb input [-o output] [-s spacing] [-r rcut] [-k kewald] [-f | --fractional_basis] [-b | --output_binary]`
* `input`: Path to input material structure in CIF format
* `output`: Path to output tabulated potential file. Default: `zeolite.ztbx`
* `spacing`: Grid spacing of the tabulate potential in angstroms.
* `rcut`: Cutoff distance of non-bonding interactions. A larger rcut results in the unit cell being duplicated more times, thus also increasing the cost and accuracy of Ewald summation.
* `kewald`: Ewald summation parameters. Default value is 3.2/rcut, use `-k 0` to turn off Ewald summation.
* `f`: Use grid direction along fractional coordinate basis, has no effect if the material has a orthorombic/cubic unit cell.
* `b`: Output the tabulated potential in binary format to be read by MCCCS-MN. If this option is turned on the `f` option is automatically applied.

## Output file format (without `-b` option)
Line 1: Unit cell information, Number of grid points along each dimension (*N_a*, *N_b*, *N_c*) and cell angles in degree (*alpha*, *beta*, *gamma*).

Line 2: Potential tabulation parameters: `spacing`, `rcut`, `kewald` and `--fractional_basis`.

Line 3: Atom types included in the unit cell.

Line 4 and after: values in tabulated potentials. Each line represents the linearized potential grid for an atom type and be reshaped into (*N_a*, *N_b*, *N_c*, 3).
