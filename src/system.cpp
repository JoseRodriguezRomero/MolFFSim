#define DEFAULT_IS_PERIODIC         {false, false, false}
#define DEFAULT_PERIODIC_SIZES      {50, 50, 50}
#define DEFAULT_SYSTEM_CHARGE       0

#define MAX_MOLEC_LABEL_SIZE        256
#define MAX_PRINT_BUFFER_SIZE       256

#define BOHR_TO_ANGSTROM            0.529177249

#include "system.hpp"

unsigned AtomicNumberFromName(const std::string &element) {
    if (element == "H") {
        return 1;
    }
    else if (element == "He") {
        return 2;
    }
    else if (element == "Li") {
        return 3;
    }
    else if (element == "Be") {
        return 4;
    }
    else if (element == "B") {
        return 5;
    }
    else if (element == "C") {
        return 6;
    }
    else if (element == "N") {
        return 7;
    }
    else if (element == "O") {
        return 8;
    }
    else if (element == "F") {
        return 9;
    }
    else if (element == "Ne") {
        return 10;
    }
    else if (element == "Na") {
        return 11;
    }
    else if (element == "Mg") {
        return 12;
    }
    else if (element == "Al") {
        return 13;
    }
    else if (element == "Si") {
        return 14;
    }
    else if (element == "P") {
        return 15;
    }
    else if (element == "S") {
        return 16;
    }
    else if (element == "Cl") {
        return 17;
    }
    else if (element == "Ar") {
        return 18;
    }
    
    return 0;
}

using namespace MolFFSim;

template<typename T>
System<T>::System() {
    is_periodic = DEFAULT_IS_PERIODIC;
    box_side_len = DEFAULT_PERIODIC_SIZES;
    system_sum_charges = 0;
}

template<typename T>
System<T>::~System() {
}

template<typename T>
void System<T>::PolarizeMolecules() {
}

template<typename T>
void System<T>::ReadInputFile(std::ifstream &input_file) {
    std::string line;
    while (std::getline(input_file, line)) {
        if (line.find("system_charge") != std::string::npos) {
            char dummy[64];
            unsigned system_charge;
            int valid_args = sscanf(line.c_str(),"%s %ud", dummy,
                                    &system_charge);
            
            if (valid_args != 2) {
                std::cout << "Incorrect arguments for \"system_charge\".";
                std::cout << std::endl;
                exit(1);
            }
            
            setSystemCharge(system_charge);
        }
        else if (line.find("is_periodic") != std::string::npos) {
            char dummy[64];
            char arg1[16], arg2[16], arg3[16];
            char *args[3] = {arg1, arg2, arg3};
            int valid_args = sscanf(line.c_str(),"%s %s %s %s", dummy, 
                                    arg1, arg2, arg3);
            
            if (valid_args != 4) {
                std::cout << "Incorrect arguments for \"is_periodic\".";
                std::cout << std::endl;
                exit(1);
            }
            
            for (unsigned i = 0; i < 3; i++) {
                if (!strcmp(args[i],"true")) {
                    is_periodic[i] = true;
                }
                else if (!strcmp(args[i],"false")) {
                    is_periodic[i] = false;
                }
                else {
                    std::cout << "Incorrect arguments for \"is_periodic\".";
                    std::cout << std::endl;
                    exit(1);
                }
            }
            
            setPeriodic(is_periodic);
        }
        else if (line.find("periodic_box_sizes") != std::string::npos) {
            char dummy[64];
            double box_sizes[3];
            int valid_args = sscanf(line.c_str(),"%s %lf %lf %lf", dummy,
                                    box_sizes, box_sizes+1, box_sizes+2);
            
            if (valid_args != 4) {
                std::cout << "Incorrect arguments for \"periodic_box_sizes\".";
                std::cout << std::endl;
                exit(1);
            }
            
            for (unsigned i = 0; i < 3; i++) {
                if (box_sizes[i] > 0) {
                    box_side_len[i] = box_sizes[i];
                }
                else {
                    std::cout << "Incorrect arguments for ";
                    std::cout << "\"periodic_box_sizes\"" << std::endl;
                    exit(1);
                }
            }
            
            setPeriodicBoxSizes(box_side_len);
        }
        else if (line.find("atomic_basis") != std::string::npos) {
            char dummy[64];
            char atom_basis_funcs_collection[256];
            
            int valid_args = sscanf(line.c_str(),"%s %s",dummy,
                                    atom_basis_funcs_collection);
            
            if (valid_args != 2) {
                std::cout << "Incorrect arguments for \"atomic_basis\".";
                std::cout << std::endl;
                exit(1);
            }
            
            this->atom_basis_funcs_collection = atom_basis_funcs_collection;
        }
        else if (line.find("xc_coefficients") != std::string::npos) {
            char dummy[64];
            char xc_coeff_collection[256];
            
            int valid_args = sscanf(line.c_str(),"%s %s",dummy,
                                    xc_coeff_collection);
            
            if (valid_args != 2) {
                std::cout << "Incorrect arguments for \"xc_coefficients\".";
                std::cout << std::endl;
                exit(1);
            }
            
            this->xc_coeff_collection = xc_coeff_collection;
        }
        else if (line == "BEGIN MOLECULE_LIST") {
            while (true) {
                if (!std::getline(input_file, line)) {
                    std::cout << "\"END MOLECULE_LIST\" is missing.";
                    std::cout << std::endl;
                    exit(1);
                }
                
                if (line.find("BEGIN MOLECULE") != std::string::npos) {
                    char molec_label[MAX_MOLEC_LABEL_SIZE];
                    char dummy1[64], dummy2[64];
                    int valid_args = sscanf(line.c_str(),"%s %s %s",
                                            dummy1,dummy2,molec_label);
                    
                    if (valid_args != 3) {
                        std::cout << "Molecule label is missing." << std::endl;
                        exit(1);
                    }
                    
                    molecule_list[molec_label] = MolFFSim::Molecule<T>();
                    
                    while (true) {
                        if (!std::getline(input_file, line)) {
                            exit(1);
                        }
                        
                        if (line == "END MOLECULE") {
                            break;
                        }
                        
                        char element_label[5];
                        double atom_x, atom_y, atom_z;
                        
                        int valid_args = 
                            sscanf(line.c_str(),"%s %lf %lf %lf",element_label,
                                   &atom_x,&atom_y,&atom_z);
                        
                        if (valid_args != 4) {
                            std::cout << "Incorrect atom in Molecule \"";
                            std::cout << molec_label << "\"." << std::endl;
                            exit(1);
                        }
                        
                        unsigned atomic_number =
                            AtomicNumberFromName(element_label);
                        
                        if (!atomic_number) {
                            std::cout << "Unsupported element in molecule ";
                            std::cout << "definition \"" << molec_label;
                            std::cout <<  "\"." << std::endl;
                            exit(1);
                        }
                        
                        atom_x /= BOHR_TO_ANGSTROM;
                        atom_y /= BOHR_TO_ANGSTROM;
                        atom_z /= BOHR_TO_ANGSTROM;
                        
                        MolFFSim::Atom<double> new_atom;
                        new_atom.setAtomicNumber(atomic_number);
                        new_atom.setPos({atom_x,atom_y,atom_z});
                        molecule_list[molec_label].addAtom(new_atom);
                    }
                }
                
                if (line == "END MOLECULE_LIST") {
                    break;
                }
            }
        }
        else if (line == "BEGIN MOLECULAR_SYSTEM") {
            while (true) {
                if (!std::getline(input_file, line)) {
                    std::cout << "\"END MOLECULAR_SYSTEM\" is missing.";
                    std::cout << std::endl;
                    exit(1);
                }
                
                if (line == "END MOLECULAR_SYSTEM") {
                    break;
                }
                
                char molec_label[MAX_MOLEC_LABEL_SIZE];
                char print_opt[32];
                char rot_opt[32];
                double dx, dy, dz;
                
                if (line.find("EULER_XYZ") != std::string::npos) {
                    double thx, thy, thz;
                    
                    int valid_args = 
                        sscanf(line.c_str(), "%s %s %s %lf %lf %lf %lf %lf %lf",
                               molec_label, print_opt, rot_opt, &dx, &dy, &dz,
                               &thx, &thy, &thz);
                    
                    if (valid_args != 9) {
                        std::cout << "Incorrect molecule specification in ";
                        std::cout << "\"MOLECULAR_SYSTEM\"." << std::endl;
                        exit(1);
                    }
                    
                    dx /= BOHR_TO_ANGSTROM;
                    dy /= BOHR_TO_ANGSTROM;
                    dz /= BOHR_TO_ANGSTROM;
                    
                    auto molec_it = molecule_list.find(molec_label);
                    if (molec_it == molecule_list.end()) {
                        std::cout << "Unknown molecular species \"";
                        std::cout << molec_label << "\"." << std::endl;
                        exit(1);
                    }
                    
                    if (!strcmp(print_opt,"PRINT")) {
                        molecules_to_print.push_back(true);
                    }
                    else if (!strcmp(print_opt,"NO_PRINT")) {
                        molecules_to_print.push_back(false);
                    }
                    else {
                        std::cout << "Print option \"" << print_opt << "\" ";
                        std::cout << "in \"MOLECULAR_SYSTEM\" is ";
                        std::cout << "incorrect." << std::endl;
                        exit(1);
                    }
                    
                    molecules.push_back(molecule_list[molec_label]);
                    molecules.back().setAnglesXYZ(thx,thy,thz);
                    molecules.back().setPos({dx,dy,dz});
                    
                    molecules.back().applyRotationAndTranslation();
                }
                else if (line.find("QUATERNION") != std::string::npos) {
                    double qw, qx, qy, qz;
                    
                    int valid_args =
                        sscanf(line.c_str(), 
                               "%s %s %s %lf %lf %lf %lf %lf %lf %lf",
                               molec_label, print_opt, rot_opt, &dx, &dy, &dz,
                               &qw, &qx, &qy, &qz);
                    
                    if (valid_args != 10) {
                        std::cout << "Incorrect molecule specification in ";
                        std::cout << "\"MOLECULAR_SYSTEM\"." << std::endl;
                        exit(1);
                    }
                    
                    dx /= BOHR_TO_ANGSTROM;
                    dy /= BOHR_TO_ANGSTROM;
                    dz /= BOHR_TO_ANGSTROM;
                    
                    auto molec_it = molecule_list.find(molec_label);
                    if (molec_it == molecule_list.end()) {
                        std::cout << "Unknown molecular species \"";
                        std::cout << molec_label << "\"." << std::endl;
                        exit(1);
                    }
                    
                    if (!strcmp(print_opt,"PRINT")) {
                        molecules_to_print.push_back(true);
                    }
                    else if (!strcmp(print_opt,"NO_PRINT")) {
                        molecules_to_print.push_back(false);
                    }
                    else {
                        std::cout << "Print option \"" << print_opt << "\" ";
                        std::cout << "in \"MOLECULAR_SYSTEM\" is ";
                        std::cout << "incorrect." << std::endl;
                        exit(1);
                    }
                    
                    molecules.push_back(molecule_list[molec_label]);
                    molecules.back().setRotQ({qw,qx,qy,qz});
                    molecules.back().setPos({dx,dy,dz});
                    
                    molecules.back().applyRotationAndTranslation();
                }
                else {
                    std::cout << "Rotation representation is incorrect or is ";
                    std::cout << "missing." << std::endl;
                    exit(1);
                }
            }
        }
    }
}

template<typename T>
void System<T>::molecGeomOptimization() {
    
}

template<typename T>
void System<T>::pointEnergyCalculation() {
    
}

template<typename T>
void System<T>::systemGeomOptimization() {
    
}

template<typename T>
void System<T>::setSystemCharge(const unsigned charge) {
    system_sum_charges = charge;
    
    auto it1_end = molecules.end();
    auto it1_begin = molecules.begin();
    for (auto it1 = it1_begin; it1 != it1_end; it1++) {
        auto it2_end = it1->Atoms().end();
        auto it2_begin = it1->Atoms().begin();
        for (auto it2 = it2_begin; it2 != it2_end; it2++) {
            system_sum_charges -= it2->AtomicNumber();
        }
    }
}

template<typename T>
void System<T>::setPeriodic(const std::vector<bool> &is_periodic) {
    this->is_periodic = is_periodic;
    
    for (auto it = molecules.begin(); it != molecules.end(); it++) {
        it->setPeriodic(is_periodic);
    }
}

template<typename T>
void System<T>::setPeriodicBoxSizes(const std::vector<double> box_side_len) {
    this->box_side_len = box_side_len;
    
    for (auto it = molecules.begin(); it != molecules.end(); it++) {
        it->setPeriodicBoxSizes(box_side_len);
    }
}

template <typename T>
std::ostream& operator<<(std::ostream &os, const MolFFSim::System<T> &system) {
    char buffer[MAX_PRINT_BUFFER_SIZE];
    for (int i = 0; i < 116; i++) {
        os << "-";
    }
    os << std::endl;
    
    snprintf(buffer, MAX_PRINT_BUFFER_SIZE,"%12s %25s %25s %25s %25s",
             "Element","x [Angstrom]","y [Angstrom]","z [Angstrom]",
             "Partial charge");
    os << buffer << std::endl;
    
    for (int i = 0; i < 116; i++) {
        os << "-";
    }
    os << std::endl;
    
    for (unsigned i = 0; i < system.Molecules().size(); i++) {
        if (!system.MoleculesToPrint()[i]) {
            continue;
        }
        
        os << system.Molecules()[i];
        
        for (int i = 0; i < 116; i++) {
            os << "-";
        }
        os << std::endl;
    }
    
    os << std::endl;
    
    snprintf(buffer,MAX_PRINT_BUFFER_SIZE, "%-20s %-15.5E [kJ/mol]",
             "Total Energy:",0.0);
    os << buffer << std::endl;
    
    snprintf(buffer,MAX_PRINT_BUFFER_SIZE, "%-20s %-15.5E [kJ/mol]",
             "Interaction Energy:",0.0);
    os << buffer << std::endl << std::endl;
    
    return os;
}

// Explicit instantiation of all the types for the class System.
template class MolFFSim::System<double>;
template class MolFFSim::System<autodiff::var>;
template class MolFFSim::System<autodiff::dual>;

template std::ostream& operator<<(std::ostream &os, const MolFFSim::System<double> &system);
template std::ostream& operator<<(std::ostream &os, const MolFFSim::System<autodiff::var> &system);
template std::ostream& operator<<(std::ostream &os, const MolFFSim::System<autodiff::dual> &system);
