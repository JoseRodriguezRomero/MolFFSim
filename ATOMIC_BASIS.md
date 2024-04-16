# Atom centered basis functions
The program approximates atomic electronic ground state eigendensities using s-type Gaussian basis functions centered at the atom's position. Mathematically, this is expressed as:
$$
\rho_a \left( r \right) = \sum_k c_{a,k} \left( \frac{\lambda_{a,k}}{\pi} \right)^{3/2} \exp \left( - \lambda_{a,k} \left \| r - r_a \right \|^2 \right),
$$
here, $a$ labels an element in the periodic table, $r_a \in \mathbb{R}^3$ denotes the vector indicating the position of the atom's nuclei, and $c_{a,k}$ and $\lambda_{a,k}$ are scalar coefficients obtained from least-squares fitting against *ab initio* calculated electronic ground state eigendensities.

For eigendensities calculated with an [Effective Core Potential (ECP)](https://manual.q-chem.com/5.3/Ch8.S10.html) basis-set, the program uses three s-type Gaussian basis functions to approximate it. Conversely, for eigendensities calculated with a full electron basis-set, nine s-type Gaussian basis functions are utilized, except for hydrogen and helium atoms, where, similar to ECP densities, only three basis functions are used. In the latter case, the last three terms in the summation account for contributions from the valence electron atomic shells.

 * `ecp_collection_1`
   * **Reference:** [*J. Phys. Chem. A,* 128, 6, 1163–1172, (2024)](https://doi.org/10.1021/acs.jpca.3c06724)
   * **Type:** ECP
   * **Coverage:** All the elements in the first three rows of the periodic table.

 * `fulle_collection_1`
   * **Reference:** [*J. Chem. Phys.* In review (2024)](https://doi.org/10.26434/chemrxiv-2024-t5tfh-v2)
   * **Type:** Full electron
   * **Coverage:** All the elements in the first three rows of the periodic table.