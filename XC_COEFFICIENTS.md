# Exchange and Correlation (XC) coefficients
To calculate the interaction energy of a system with more than one atom, we first compute the interaction energy as if the entire system's wave function were a [Hartree product](http://vergil.chemistry.gatech.edu/notes/hf-intro/node3.html) of each atom's wave function, treating them as if they were non-interacting. This initial calculation yields what we refer to as the "naïve model."

To address the inaccuracies stemming from this oversimplified model, we introduce a set of fictitious electronic clouds. These clouds enable us to compute an exchange and correlation component, designed solely to correct for these inaccuracies. Each cloud is fully parametrized based on the non-interacting ground state eigendensity of each atom in the system, achieved through successive applications of differential operators.

Moreover, the exchange and correlation component is decomposed into four parts. Two parts concern electron-electron and electron-nuclei interactions with spherical symmetry in the fictitious clouds, while the other two involve cylindrical symmetry along the axis between pairs of atoms. For a system with only two atoms, labeled as $a$ and $b$, this exchange and correlation component can be expressed as:
$$
\text{XC}_{ab} = \sum_{k = 1}^K \left( - 1 \right)^k \left( A \frac{\text{XC}_{ab}^{\text{EE} \left( k \right)}}{\left( 2 k \right) !}  + B \frac{\text{XC}_{ab}^{\text{EN} \left( k \right)}}{\left( 2 k \right) !} - C \frac{\text{XC}_{ab}^{\text{ED} \left( k \right)}}{\left( 2 k + 1 \right) !} - D \frac{\text{XC}_{ab}^{\text{ND} \left( k \right)}}{\left( 2 k + 1 \right) !} \right),
$$
here, $K$ denotes the truncation order of the series (set to ten in this context), and $A$, $B$, $C$, and $D$ represent real-valued coefficients fitted against high-accuracy *ab initio* reference calculations, specifically using the CCSD(T)/cc-pVTZ level of theory. The respective functionals are defined as:
$$
\begin{align}
\text{XC}_{ab}^{\text{EE} \left( k \right)} &= \frac{1}{2} \int_{\mathbb{R}^3} \int_{\mathbb{R}^3} \left( \frac{\rho_a \left( r \right) \nabla^{2 k}_{r'} \rho_b \left( r' \right)}{\left \lVert r - r' \right \rVert} + \frac{\rho_b \left( r \right) \nabla^{2 k}_{r'} \rho_a \left( r' \right)}{\left \lVert r - r' \right \rVert}\right) \mathrm{d} r \mathrm{d} r' \\
\text{XC}_{ab}^{\text{ED} \left( k \right)} &= \frac{1}{2} \int_{\mathbb{R}^3} \int_{\mathbb{R}^3} \left( \frac{\rho_a \left( r \right) r_{ba}}{\left \lVert r_{ba} \right \rVert} \cdot \frac{\nabla^{2 k + 1}_{r'} \rho_b \left( r' \right)}{\left \lVert r - r' \right \rVert} + \frac{\rho_b \left( r \right) r_{ab}}{\left \lVert r_{ab} \right \rVert} \cdot \frac{\nabla^{2 k + 1}_{r'} \rho_a\left( r' \right)}{\left \lVert r - r' \right \rVert} \right) \mathrm{d} r \mathrm{d} r' \\
\text{XC}_{ab}^{\text{EN} \left( k \right)} &= \int_{\mathbb{R}^3} \left( \frac{Z_a \nabla^{2 k} \rho_b \left( r \right)}{\left \lVert r - r_a \right \rVert} +  \frac{Z_b \nabla^{2 k} \rho_a \left( r \right)}{\left \lVert r - r_b \right \rVert} \right) \mathrm{d} r \\
\text{XC}_{ab}^{\text{ND} \left( k \right)} &= \int_{\mathbb{R}^3} \left( \frac{ Z_a r_{ba}}{\left \lVert r_{ba} \right \rVert} \cdot \frac{\nabla^{2 k + 1} \rho_b \left( r \right)}{\left \lVert r - r_a \right \rVert} + \frac{ Z_b r_{ab}}{\left \lVert r_{ab} \right \rVert} \cdot \frac{\nabla^{2 k + 1} \rho_a \left( r \right)}{\left \lVert r - r_b \right \rVert} \right) \mathrm{d} r,
\end{align}
$$
where $\rho_a, \rho_b \in L^1 \left( \mathbb{R^3} \right)$ are the non-interacting ground state eigendensities of the atoms, $Z_a$ and $Z_b$ denote their atomic numbers, $r_a$ and $r_b$ represent the position vectors of the atoms, and $r_{ab} = r_b - r_a = -  r_{ba}$ represents the distance vectors between the two atoms.

To accommodate partial charge transfer between the atoms, the valence electron contributions of the non-interacting ground state eigendensities of each atom are scaled by a polarization scalar. These scalars are determined by minimizing a function identical in form to the energy calculation formula, subject to the constraint that the integral of the densities remains equal to the total number of electrons in the entire system. Additionally, the exchange and correlation coefficients used in the minimization problem differ from those used to compute the energies, resulting in a total of eight fitted parameters to account for exchange and correlation.

 * `ecp_collection_1`
   * **Reference:** [*J. Chem. Phys.* In review (2024)](https://doi.org/10.26434/chemrxiv-2024-t5tfh-v2)
   * **Fitted with:** The atomic basis functions in `ecp_collection_1`.

 * `fulle_collection_1`
   * **Reference:** [*J. Chem. Phys.* In review (2024)](https://doi.org/10.26434/chemrxiv-2024-t5tfh-v2)
   * **Fitted with:** The atomic basis functions in `fulle_collection_1`.