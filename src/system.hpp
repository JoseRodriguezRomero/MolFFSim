#ifndef CLASS_MOLECULES_SYSTEM
#define CLASS_MOLECULES_SYSTEM

#include <stdio.h>
#include <lbfgs.h>

#include <cmath>
#include <thread>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <functional>
#include <unordered_map>

#include <autodiff/forward/dual.hpp>

#include "auxiliary_functions.hpp"
#include "molecule.hpp"

namespace MolFFSim {
template<typename T>
class System {
private:
    // All the molecules in the system.
    std::vector<Molecule<T>> molecules;
    std::vector<bool> molecules_to_print;
    std::vector<std::string> names_molecules;
    
    // This is are auxiliary vectors, and should never be return in any of the
    // interfaces of this class.
    std::vector<Atom<T>*> atoms_molecules;
    
    // All the type of molecules in the system. The keys are the user-assigned
    // ID's (in the form of strings) given in the input file.
    std::unordered_map<std::string,MolFFSim::Molecule<T>> molecule_list;
    std::unordered_map<std::string,unsigned> molecule_instances;
    std::unordered_map<std::string,double> monomer_energies;
    
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
    
    Eigen::Vector<T,Eigen::Dynamic> sys_params;
    const std::thread::id t_id = std::this_thread::get_id();
    
public:
    System();
    ~System();
    
    const std::vector<Molecule<T>>& Molecules() const { return molecules; }
    const std::vector<bool>& MoleculesToPrint() const {
        return molecules_to_print;
    }
    
    void ReadInputFile(std::ifstream &input_file);
    
    // A vector with 8*N entries, where N is the number of molecules in the
    // system. The first three entries are the displacement vector, the
    // next four are the components of the rotaiton quaternion. The last
    // entry per molecule is a Lagrange multiplier to make the quaternions
    // normal and avoif singularies.
    inline Eigen::Vector<T,Eigen::Dynamic> SysParams() const {
        return sys_params;
    }

    void SetSysParams(const Eigen::Vector<T,Eigen::Dynamic> &sys_params);
    
    inline T
    EnergyFromParams(const Eigen::Vector<T,Eigen::Dynamic> &sys_params) {
        SetSysParams(sys_params);
        
        T energy = SystemEnergy();
        for (unsigned i = 0; i < molecules.size(); i++) {
            T aux_sum = -1.0;
            aux_sum += pow(sys_params(i*7 + 3),2);
            aux_sum += pow(sys_params(i*7 + 4),2);
            aux_sum += pow(sys_params(i*7 + 5),2);
            aux_sum += pow(sys_params(i*7 + 6),2);
            energy += pow(aux_sum,2);
        }
        
        return energy;
    }
    
    Eigen::Vector<T,Eigen::Dynamic> 
    GradEnergyFromParams(const Eigen::Vector<T,Eigen::Dynamic> &sys_params);
        
    // Use this for system-wide a geometry optimization, and point-energy
    // calculations of the system as a whole.
    void PolarizeMolecules();

    T SystemEnergy();
    T SystemInteractionEnergy();
    
    const std::vector<std::string> ListMoleculeTypes() const;
    const unsigned MoleculeInstances(const std::string &molec_name) const;
    const double MoleculeMonomerEnergy(const std::string &molec_name) const;
    
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
    
    int OptimizeGeometry(std::ostream &os = std::cout);
};

}

template <typename T>
std::ostream& operator<<(std::ostream &os, const MolFFSim::System<T> &system);

#endif
