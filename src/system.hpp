#ifndef CLASS_MOLECULES_SYSTEM
#define CLASS_MOLECULES_SYSTEM

#include <stdio.h>

#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <unordered_map>

#include <autodiff/reverse/var.hpp>
#include <autodiff/forward/dual.hpp>

#include "auxiliary_functions.hpp"
#include "molecule.hpp"

namespace MolFFSim {
template<typename T>
class System {
private:
    // All the molecules in the system.
    std::vector<Molecule<T>> molecules;
    std::vector<double> monomer_energies;
    std::vector<bool> molecules_to_print;
    std::vector<std::string> names_molecules;
    
    // All the type of molecules in the system. The keys are the user-assigned
    // ID's (in the form of strings) given in the input file.
    std::unordered_map<std::string,MolFFSim::Molecule<T>> molecule_list;
    
    // Basis-function parametrization for the non-interacting ground states
    // eigendensities of all the elements of interest. The atomic number of
    // each element are used as hash keys.
    std::unordered_map<unsigned,std::vector<double>> elem_c_coeff;
    std::unordered_map<unsigned,std::vector<double>> elem_lambda_coeff;
    bool is_ecp;
    
    std::unordered_map<unsigned,std::vector<double>> elem_xc_coeff;
    MolFFSim::XCRules xc_rule = rule1;
    
    std::string xc_coeff_collection;            // Name of the collection.
    std::string atom_basis_collection;          // Name of the collection.
    
    std::vector<bool> is_periodic;
    std::vector<double> box_side_len;
    
    unsigned system_sum_charges;
    
public:
    System();
    ~System();
    
    const std::vector<Molecule<T>>& Molecules() const { return molecules; }
    const std::vector<bool>& MoleculesToPrint() const {
        return molecules_to_print;
    }
    
    void ReadInputFile(std::ifstream &input_file);
    
    // Use this for system-wide a geometry optimization, and point-energy
    // calculations of the system as a whole.
    void PolarizeMolecules();

    T SystemEnergy();
    T SystemInteractionEnergy();
    T SystemEnergyFromGeom(const std::vector<T> angles_and);
    
    // Use this to optimize the isolated (monomer) geometry of the molecules.
    void MonomerPolarizeMolecules();
    void getMonomerEnergies();
    
    void setSystemCharge(const unsigned charge);
    void setPeriodic(const std::vector<bool> &is_periodic);
    void setPeriodicBoxSizes(const std::vector<double> box_side_len);
    
    void readXCCoefficients(const std::string &xc_coeff_collection);
    void readAtomicBasisSet(const std::string &atom_basis_collection);
    
    void setECP();
    void setFullE();
    void setXCRule();
};

}

template <typename T>
std::ostream& operator<<(std::ostream &os, const MolFFSim::System<T> &system);

#endif
