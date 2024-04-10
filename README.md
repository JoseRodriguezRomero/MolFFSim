![Static Badge](https://img.shields.io/badge/Linux-untested-yellow?logo=linux&logoColor=white)
![Static Badge](https://img.shields.io/badge/MacOS-runs-green?logo=apple&logoColor=white) 
![Static Badge](https://img.shields.io/badge/Windows-unsupported-red?logo=windows&logoColor=white) \
![Static Badge](https://img.shields.io/badge/g%2B%2B-compiles-green?logo=gnu&logoColor=white)
![Static Badge](https://img.shields.io/badge/clang-compiles-green?logo=llvm&logoColor=white)
![Static Badge](https://img.shields.io/badge/msvc%2B%2B-untested-yellow?logo=microsoft&logoColor=white) \
![Static Badge](https://img.shields.io/badge/OpenMP-supported-green)
![Static Badge](https://img.shields.io/badge/MPI-unsupported-red) \
![Static Badge](https://img.shields.io/badge/license-GPL--3.0-blue)
![Static Badge](https://img.shields.io/badge/DOI-10.26434%2Fchemrxiv--2024--t5tfh--v2-blue?link=https%3A%2F%2Fchemrxiv.org%2Fengage%2Fchemrxiv%2Farticle-details%2F66033bf59138d23161602c69)
# MolFFSim
MolFFSim is an open-source program coded in C++ specifically designed for performing point-energy calculations of systems with one or many molecules, as well as optimizing their geometry by minimizing the energy. In the latter, all molecules are treated as rigid objects. The theoretical framework for the calculations in this program is rooted in the Orbital Free Density Functional Theory (OFDFT) approach, as detailed in the following references:
  * [*J. Phys. Chem. A,* 128, 6, 1163–1172 (2024)](https://pubs.acs.org/doi/10.1021/acs.jpca.3c06724)
  * [*J. Chem. Phys.* In review (2024)](https://doi.org/10.26434/chemrxiv-2024-t5tfh-v2)

Parallelization support is exclusively achieved through OpenMP. Exact values of gradients and Hessians of interaction energies are obtained using algorithmic differentiation via [autodiff](https://autodiff.github.io/). Both forward and reverse modes are fully supported, without the need of recompiling the source code. All the linear algebra necessary for calculations within this program is handled using [Eigen](https://gitlab.com/libeigen/eigen).

## Compiling and installing
In principle, any C++ compiler that supports C++17 can compile this program. However, thus far, it has been tested only with GCC and Clang compilers. Before proceeding with the compilation and installation, ensure that CMake is installed, along with the C++ libraries [autodiff](https://autodiff.github.io/) and [Eigen](https://gitlab.com/libeigen/eigen). Detailed instructions on compiling and installing these libraries are available on their respective GitHub/GitLab repositories.

Once these prerequisites are fulfilled, follow these steps to compile and install the program:
```Bash
cd Downloads
git clone https://github.com/JoseRodriguezRomero/MolFFSim.git
cd MolFFSim
mkdir build && cd build
cmake ..
make && sudo make install
```

## Running a simulation
After installing the program on your computer, conducting a point-energy calculation for a system is straightforward. Simply execute the command `MolFFSim <input_file>` in a Bash or Z shell.

### Input files
To run this program without errors at least the name of an existing input file has to be provided. The name of the input file must _always_ be provided as the first argument. Furthermore, the content of this input file must be formatted according to the rules described in this section. For clarity, two example of valid input files are provided with the source code. 

The structure of an input files can be separated in three main parts, that we will now describe.

#### 1. System Configuration
 * `system_charge` Specifies the overall charge of the system.
 * `is_periodic` Indicates whether the system is periodic or not, with three boolean values representing periodicity along the *x*, *y*, and *z* axes respectively.
 * `periodic_box_sizes` Specifies the length of the sides of the periodic box (in Ångström), centered at the origin, along the *x*, *y*, and *z* axes.
 * `atomic_basis` Identifies the collection of basis functions used for atomic representations.
 * `xc_coefficients` Identifies the collection of exchange-correlation coefficients used.

#### 2. Molecule List
The molecule list contains descriptions of individual molecules in the system. Each molecule is defined by its name followed by its atomic composition and coordinates.

 * `BEGIN MOLECULE_LIST` Marks the beginning of the molecule list section.
 * `BEGIN MOLECULE <molecule_name>` Marks the beginning of a molecule with its name.
 * Atomic symbols followed by their *x*, *y*, and *z* coordinates (in Ångström) define each atom within the molecule.
 * `END MOLECULE` Marks the end of a molecule's description.
 * `END MOLECULE_LIST` Marks the end of the molecule list section.

#### 3. Molecular System Configuration
This part defines the configuration of the molecular system, including how individual molecules are positioned and whether they are printed or not.

 * `BEGIN MOLECULAR_SYSTEM` Marks the beginning of the molecular system configuration.
 * Each line describes a molecule within the system, including its name, print status, rotation type, and orientation parameters. For Euler angles the last three arguments describe the rotation angles, and for quaternions the last four arguments describe its *w*, *x*, *y* and *z* real and complex components. The first three numerical arguments correspond to the distance vector (in Ångström) at which the molecule is displaced.
 * `END MOLECULAR_SYSTEM` Marks the end of the molecular system configuration.

### Point-energy calculations
By default, all calculations performed are point-energy calculations. This means that when invoking the `MolFFSim` executable in a console with just the input file's name, the calculation type is automatically assumed to be a point-energy calculation. This is equivalent to executing the command:
```Bash
MolFFSim <input_file> point_energy
```

The output of this type of calculation provides both the total and interaction energies of the entire system, along with the coordinates and partial charges of the molecules. However, molecules marked with `NO_PRINT` in the `MOLECULAR_SYSTEM` section of the input file will be excluded from the output.

### Geometry optimization of a molecule
> [!CAUTION]
> This type of calculation is experimental, and with any of the predefined sets of exchange and correlation coefficients, the results obtained lack physical relevance.

The second type of calculation supported in this program involves optimizing the intramolecular coordinates of a given molecule. For this type of calculation, the entire `MOLECULAR_SYSTEM` section of the input file can be omitted, and if present, the program simply ignores it. To run this type of calculation, use the following command: 
```Bash
MolFFSim <input_file> molec_geom_optim
```
In the output of this calculation, the Cartesian coordinates and partial charges of every atom of the molecule are printed. If more than one molecule is specified in the `MOLECULE_LIST` section of the input file, the optimized geometry of every molecule in the list is printed. The calculation of each molecule is independent of the others, making it [embarrassingly parallelizable](https://en.wikipedia.org/wiki/Embarrassingly_parallel#:~:text=In%20parallel%20computing%2C%20an%20embarrassingly,a%20number%20of%20parallel%20tasks.).

### Geometry optimization of a system with two or more molecules
This program also supports system-wide geometry optimization for ensembles containing two or more molecules. In this mode, each molecule in the system is treated as a rigid object, with translations and rotations applied uniformly to their atoms. To initiate this calculation, use the following command:
```Bash
MolFFSim <input_file> system_geom_optim
```
The output format of this calculation mirrors that of a point-energy calculation. 

### Forward and reverse mode automatic differentiation
To model the [Mulliken charges](https://chem.libretexts.org/Bookshelves/Physical_and_Theoretical_Chemistry_Textbook_Maps/Book%3A_Quantum_States_of_Atoms_and_Molecules_(Zielinksi_et_al)/10%3A_Theories_of_Electronic_Molecular_Structure/10.07%3A_Mulliken_Populations) of a molecular system, a fictitious self-energy functional is minimized. This functional mirrors the form of the functional used for modeling the system energy. Minimizing this problem simplifies to minimizing a quadratic function, thus solving a system of linear equations. The size of this system scales linearly with the number of atoms in the system.

To derive expressions for the gradient and Hessian, knowledge of the matrix inverse associated with the aforementioned minimization problem is required. Therefore, algorithmic differentiation is used.

To opt for forward mode automatic differentiation during a geometry optimization of a system, use the following command:
```Bash
MolFFSim <input_file> system_geom_optim forward_mode
```

For reverse mode, use:
```Bash
MolFFSim <input_file> system_geom_optim reverse_mode
```

If not specified, reverse mode automatic differentiation is the default. Additionally, for point-energy calculations, gradients or Hessians are unnecessary. Thus, the choice between `forward_mode` or `reverse_mode` is inconsequential.

## Atom centered basis functions
The program's core functionality revolves around determining energies and partial charges. These are computed based on the non-interacting ground state eigendensities of each atom in the system. These densities are approximated by a linear combination of s-type Gaussian basis functions centered at the atomic coordinates.

The electronic densities from *ab initio* calculations utilizing an [Effective Core Potential (ECP)](https://manual.q-chem.com/5.3/Ch8.S10.html) basis-set are approximated with three s-type Gaussian basis functions at each atomic site. In contrast, for electronic densities obtained from full electron basis-set calculations, nine basis functions are utilized, except for hydrogen and helium atoms, where only three basis functions suffice. The parameters governing these basis functions depend solely on the argument provided to `atomic_basis` in the input file.

For further details, please consult `ATOMIC_BASIS.md`.

## Exchange and Correlation (XC) coefficients
The energy functional in this program comprises two distinct components:
 * A naive model component, which relies on the electronic density derived from a [Hartree product](http://vergil.chemistry.gatech.edu/notes/hf-intro/node3.html) wave function. This model assumes that electrons "belonging" to different atoms behave independently, albeit contradicting the Pauli exclusion principle.
 * An exchange and correlation component, designed to rectify inaccuracies introduced by the naive model in both the potential and kinetic energy distributions of the electrons.

Given the spherical symmetry of ground state eigendensities, exact solutions exist for the energy contributions of the naive model. Moreover, fictitious electronic clouds, easily constructed via successive applications of differential operators, can be unambiguously mapped from non-interacting ground state eigendensities. Exact solutions also exist for integrals resembling Coulomb's law using these fictitious integrals.

Leveraging this understanding, the exchange and correlation functional is divided into four parts: two for electron-electron and electron-nuclei interactions, where the fictitious electron cloud exhibits spherical symmetry, and another two for electron-electron and electron-nuclei interactions, where the fictitious cloud has cylindrical symmetry. However, to fully parameterize the force field model due to different coefficients used for modeling the Mulliken charges, a total of eight fitting parameters are employed.

For further details, please consult `XC_COEFFICIENTS.md`.

## Citing this program
If you find this program useful please consider citing this repository with the following Bibtex entry:
