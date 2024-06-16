#ifndef CLASS_MOLECULES_SYSTEM
#define CLASS_MOLECULES_SYSTEM

#include <stdio.h>
#include <lbfgs.h>

#include <cmath>
#include <chrono>
#include <vector>
#include <string>
#include <utility>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <unordered_map>

#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>

#include "auxiliary_functions.hpp"
#include "molecule.hpp"

namespace MolFFSim {
template<typename T>
class System {
private:
    // All the molecules in the system.
    mutable std::vector<Molecule<T>> molecules;
    std::vector<bool> molecules_to_print;
    std::vector<std::string> names_molecules;
    
    std::vector<Molecule<T>*> relax_molecules;
    std::vector<Molecule<T>*> freeze_molecules;
    
    // This is are auxiliary vectors, and should never be return in any of the
    // interfaces of this class.
    std::vector<Atom<T>*> atoms_molecules;
    
    // unrepeated pairs, e.g. ((i,j) == (j,i)) is true.
    std::vector<std::pair<unsigned,unsigned>> atom_pairs;
    
    // All the type of molecules in the system. The keys are the user-assigned
    // ID's (in the form of strings) given in the input file.
    std::unordered_map<std::string,double> monomer_energies;
    std::unordered_map<std::string,unsigned> molecule_instances;
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
    
    unsigned system_charge;
    double uncorr_energy;
    
    unsigned n_cores;
    Eigen::Vector<T,Eigen::Dynamic> sys_params;
    
    std::ostream *output_stream;
    std::string backup_filename1;
    std::string backup_filename2;
    std::chrono::steady_clock::time_point backup_clock;
    
    mutable Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> pol_mat;
    mutable Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> aux_vec_mat;
    mutable Eigen::Vector<T,Eigen::Dynamic> energy_vec;
    
    lbfgs_parameter_t lbfgs_settings;
    
public:
    System();
    ~System();
    
    const std::vector<Molecule<T>>& Molecules() const { return molecules; }
    const std::vector<bool>& MoleculesToPrint() const {
        return molecules_to_print;
    }
    
    void ReadInputFile(std::ifstream &input_file);
    
    // A vector with 7*N entries, where N is the number of molecules in the
    // system. The first three entries are the displacement vector, the
    // next four are the components of the rotation quaternion.
    inline Eigen::Vector<T,Eigen::Dynamic> SysParams() const {
        return sys_params;
    }

    // ----------------------------------------------------------
    // --------  INTENDED FOR RIGID OBJECT MANIPULATION ---------
    // ----------------------------------------------------------
    // Every molecule is treated a rigid object here
    
    void SetSysParams(const Eigen::Vector<T,Eigen::Dynamic> &sys_params);
    
    inline T
    EnergyFromParams(const Eigen::Vector<T,Eigen::Dynamic> &sys_params) {
        SetSysParams(sys_params);
        
        T energy = SystemInteractionEnergy();
        for (unsigned i = 0; i < relax_molecules.size(); i++) {
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
    
    // ----------------------------------------------------------
    // --------  INTENDED FOR RIGID OBJECT MANIPULATION ---------
    // ----------------------------------------------------------
    
    // ----------------------------------------------------------
    // --------     INTENDED FOR CALCULATING FORCES     ---------
    // ----------------------------------------------------------
    // The whole system is manipulated on an atom-per-atom basis here
        
    void SetAtomsCoords(const Eigen::Vector<T,Eigen::Dynamic> &atom_coords);
    
    inline T
    EnergyFromAtomsCoords(const Eigen::Vector<T,
                          Eigen::Dynamic> &atom_coords) {
        SetAtomsCoords(atom_coords);
        return SystemEnergy();
    }
    
    Eigen::Vector<T,Eigen::Dynamic>
    GradEnergyFromAtomsCoords(const Eigen::Vector<T,
                              Eigen::Dynamic> &atom_coords);
    
    void printAtomForces(std::ostream &os = std::cout);
    
    // ----------------------------------------------------------
    // --------     INTENDED FOR CALCULATING FORCES     ---------
    // ----------------------------------------------------------
        
    // Use this for system-wide a geometry optimization, and point-energy
    // calculations of the system as a whole.
    void PolarizeMolecules() const;

    T SystemEnergy() const;
    inline T SystemInteractionEnergy() const {
        return SystemEnergy() - uncorr_energy;
    }
    
    std::vector<std::string> ListMoleculeTypes() const;
    unsigned MoleculeInstances(const std::string &molec_name) const;
    double MoleculeMonomerEnergy(const std::string &molec_name) const;
    
    // Use this to optimize the isolated (monomer) geometry of the molecules.
    void MonomerPolarizeMolecules() const;
    void getMonomerEnergies();
    
    void setPeriodic(const std::vector<bool> &is_periodic);
    void setPeriodicBoxSizes(const std::vector<double> box_side_len);
    
    void readXCCoefficients(const std::string &xc_coeff_collection);
    void readAtomicBasisSet(const std::string &atom_basis_collection);
    
    void setECP();
    void setFullE();
    void setXCRule();
    
    int OptimizeGeometry(std::ostream &os = std::cout);
    inline std::ostream& OStream() { return *output_stream; }
    
    int OptimizeMoleculeGeometries(std::ostream &os = std::cout);
    void printInputSettings(std::ostream &os = std::cout) const;
    
    void printBackups();
    std::string BackupFile1() const { return backup_filename1; }
    std::string BackupFile2() const { return backup_filename2; }
    
    inline void setBackupFile1(const std::string &backup_filename1) {
        this->backup_filename1 = backup_filename1;
    }
    
    inline void setBackupFile2(const std::string &backup_filename2) {
        this->backup_filename2 = backup_filename2;
    }
    
    inline std::chrono::steady_clock::time_point& backupClock(){
        return backup_clock;
    };
};

}

template <typename T>
std::ostream& operator<<(std::ostream &os, MolFFSim::System<T> &system);

#endif
