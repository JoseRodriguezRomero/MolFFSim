#define DEFAULT_IS_PERIODIC         {false, false, false}
#define DEFAULT_PERIODIC_SIZES      {50, 50, 50}
#define DEFAULT_SYSTEM_CHARGE       0

#define XC_COEFFS_BASE_DIR          "/usr/local/share/MolFFSim/XC_COEFFICIENTS/"
#define ATOMIC_BASIS_BASE_DIR       "/usr/local/share/MolFFSim/ATOMIC_BASIS/"

#define MAX_MOLEC_LABEL_SIZE        256
#define MAX_PRINT_BUFFER_SIZE       256

#define BOHR_TO_ANGSTROM            0.529177249

#include "system.hpp"

using namespace MolFFSim;

template<typename T>
System<T>::System() {
    is_periodic = DEFAULT_IS_PERIODIC;
    box_side_len = DEFAULT_PERIODIC_SIZES;
    system_sum_charges = 0;
}

template<typename T>
System<T>::~System() {
//    std::cout << std::endl;
//    for (auto it1 = molecule_list.begin(); it1 != molecule_list.end(); it1++) {
//        auto it2_end = it1->second.Atoms().cend();
//        auto it2_begin = it1->second.Atoms().cbegin();
//        for (auto it2 = it2_begin; it2 != it2_end; it2++) {
//            char buffer[256];
//            snprintf(buffer,256,"%5s %15.5lf",
//                     LabelFromAtomicNumber(it2->AtomicNumber()).c_str(),
//                     it2->CloudCCoeffs().at(0));
//            std::cout << buffer << std::endl;
//            
//            auto it3_end = it2->CloudCCoeffs().cend();
//            auto it3_begin = it2->CloudCCoeffs().cbegin() + 1;
//            for (auto it3 = it3_begin; it3 != it3_end; it3++) {
//                snprintf(buffer,256,"%5s %15.5lf", " ", *it3);
//                std::cout << buffer << std::endl;
//            }
//        }
//    }
}

template<typename T>
void System<T>::ReadInputFile(std::ifstream &input_file) {
    molecules.clear();
    monomer_energies.clear();
    molecules_to_print.clear();
    names_molecules.clear();
    molecule_instances.clear();
    atoms_molecules.clear();
    molecule_list.clear();
    elem_c_coeff.clear();
    elem_lambda_coeff.clear();
    elem_xc_coeff.clear();
    system_sum_charges = 0;
    
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
            char atom_basis_collection[256];
            
            int valid_args = sscanf(line.c_str(),"%s %s",dummy,
                                    atom_basis_collection);
            
            if (valid_args != 2) {
                std::cout << "Incorrect arguments for \"atomic_basis\".";
                std::cout << std::endl;
                exit(1);
            }
            
            this->atom_basis_collection = atom_basis_collection;
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
                            AtomicNumberFromLabel(element_label);
                        
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
                    names_molecules.push_back(molec_label);
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
                    names_molecules.push_back(molec_label);
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
    
    readXCCoefficients(xc_coeff_collection);
    readAtomicBasisSet(atom_basis_collection);
    
    if (is_ecp) {
        setECP();
    }
    else {
        setFullE();
    }
    
    setXCRule();
    
    for (auto it1 = molecules.begin(); it1 != molecules.end(); it1++) {
        auto it2_end = it1->AtomsRotAndTrans().end();
        auto it2_begin = it1->AtomsRotAndTrans().begin();
        for (auto it2 = it2_begin; it2 != it2_end; it2++) {
            atoms_molecules.push_back(&(*it2));
        }
    }
    
    getMonomerEnergies();
    for (std::string molec_name : names_molecules) {
        if (molecule_instances.find(molec_name) == molecule_instances.end()) {
            molecule_instances[molec_name] = 1;
        }
        else {
            molecule_instances[molec_name]++;
        }
    }
}

template<typename T>
void System<T>::readXCCoefficients(const std::string &xc_coeff_collection) {
    this->xc_coeff_collection = xc_coeff_collection;
    
    const std::string file_name =
        XC_COEFFS_BASE_DIR + xc_coeff_collection + ".txt";
    
    std::ifstream data_file;
    data_file.open(file_name, std::fstream::in);
    
    std::string line;
    while (std::getline(data_file, line)) {
        if (line.find("comb_rule") != std::string::npos) {
            char dummy[32], basis_type[32];
            int num_args = sscanf(line.c_str(),"%s %s",dummy,basis_type);
            
            if (num_args != 2) {
                std::cout << xc_coeff_collection << ".txt is tampered!";
                std::cout << std::endl;
                exit(1);
            }
            
            if (!strcmp(basis_type,"rule_1")) {
                xc_rule = rule1;
            }
            else if (!strcmp(basis_type,"rule_2")) {
                xc_rule = rule2;
            }
            else if (!strcmp(basis_type,"rule_3")) {
                xc_rule = rule3;
            }
            else {
                std::cout << xc_coeff_collection << ".txt is tampered!";
                std::cout << std::endl;
                exit(1);
            }
            
            continue;
        }
        
        char elem_label[8];
        double xc_coeffs[8];
        int num_args = sscanf(line.c_str(),"%s %lf %lf %lf %lf %lf %lf %lf %lf",
                              elem_label, xc_coeffs + 0, xc_coeffs + 1,
                              xc_coeffs + 2, xc_coeffs + 3, xc_coeffs + 4,
                              xc_coeffs + 5, xc_coeffs + 6, xc_coeffs + 7);
        
        if (num_args != 9) {
            std::cout << xc_coeff_collection << ".txt is tampered!";
            std::cout << std::endl;
            exit(1);
        }
        
        unsigned atomic_num = AtomicNumberFromLabel(elem_label);
        if (atomic_num == 0) {
            std::cout << atom_basis_collection << ".txt is tampered!";
            std::cout << std::endl;
            exit(1);
        }
        
        elem_xc_coeff[atomic_num] = std::vector<double>();
        elem_xc_coeff[atomic_num].insert(elem_xc_coeff[atomic_num].begin(),
                                         xc_coeffs, xc_coeffs + 8);
    }
    
    data_file.close();
}

template<typename T>
void System<T>::readAtomicBasisSet(const std::string &atom_basis_collection) {
    this->atom_basis_collection = atom_basis_collection;
    
    const std::string file_name = 
        ATOMIC_BASIS_BASE_DIR + atom_basis_collection + ".txt";
    
    std::ifstream data_file;
    data_file.open(file_name, std::fstream::in);
        
    std::string line;
    std::string current_elem = ".";
    while (std::getline(data_file, line)) {
        if (line.find("basis_type") != std::string::npos) {
            char dummy[32], basis_type[32];
            int num_args = sscanf(line.c_str(),"%s %s",dummy,basis_type);
            
            if (num_args != 2) {
                std::cout << atom_basis_collection << ".txt is tampered!";
                std::cout << std::endl;
                exit(1);
            }
            
            if (!strcmp(basis_type,"FullE")) {
                is_ecp = false;
            }
            else if (!strcmp(basis_type,"ECP")) {
                is_ecp = true;
            }
            else {
                std::cout << atom_basis_collection << ".txt is tampered!";
                std::cout << std::endl;
                exit(1);
            }
            
            continue;
        }
        
        char elem_label[8];
        double c_coeffs[3], lambda_coeffs[3];
        int num_args = sscanf(line.c_str(),"%s %lf %lf %lf %lf %lf %lf",
                              elem_label, c_coeffs,lambda_coeffs,
                              c_coeffs + 1, lambda_coeffs + 1,
                              c_coeffs + 2, lambda_coeffs + 2);
        
        if (num_args != 7) {
            std::cout << atom_basis_collection << ".txt is tampered!";
            std::cout << std::endl;
            exit(1);
        }
        
        if (strcmp(elem_label,".")) {
            current_elem = elem_label;
        }
        
        if (current_elem == ".") {
            std::cout << atom_basis_collection << ".txt is tampered!";
            std::cout << std::endl;
            exit(1);
        }
        
        std::vector<double> new_c_coeffs;
        new_c_coeffs.insert(new_c_coeffs.begin(),c_coeffs,c_coeffs+3);
        
        std::vector<double> new_lambda_coeffs;
        new_lambda_coeffs.insert(new_lambda_coeffs.begin(),lambda_coeffs,
                                 lambda_coeffs+3);
        
        unsigned atomic_num = AtomicNumberFromLabel(current_elem);
        if (atomic_num == 0) {
            std::cout << atom_basis_collection << ".txt is tampered!";
            std::cout << std::endl;
            exit(1);
        }
        
        if (elem_c_coeff.find(atomic_num) == elem_c_coeff.end()) {
            elem_c_coeff[atomic_num] = std::vector<double>();
            elem_lambda_coeff[atomic_num] = std::vector<double>();
        }
        
        for (unsigned i = 0; i < 3; i++) {
            elem_c_coeff[atomic_num].push_back(new_c_coeffs[i]);
            elem_lambda_coeff[atomic_num].push_back(new_lambda_coeffs[i]);
        }
    }
    
    data_file.close();
    
    for (auto it1 = molecules.begin(); it1 != molecules.end(); it1++) {
        it1->createElectronClouds(elem_c_coeff,elem_lambda_coeff);
        it1->setXCCoefficients(elem_xc_coeff);
    }
    
    for (auto it1 = molecule_list.begin(); it1 != molecule_list.end(); it1++) {
        it1->second.createElectronClouds(elem_c_coeff,elem_lambda_coeff);
        it1->second.setXCCoefficients(elem_xc_coeff);
    }
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
            system_sum_charges -= it2->EffAtomicNumber();
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

template<typename T>
void System<T>::PolarizeMolecules() {
    unsigned n_atoms = atoms_molecules.size();
    Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> pol_mat(n_atoms+1,n_atoms+1);
    Eigen::Matrix<T,Eigen::Dynamic,1> vec_mat(n_atoms+1);
        
    for (unsigned i = 0; i < n_atoms; i++) {
        const Atom<T> *atom_i = atoms_molecules[i];
        for (unsigned j = i + 1; j < n_atoms; j++) {
            const Atom<T> *atom_j = atoms_molecules[j];
            T mat_elem(0.0);
            T vec_elem_i(0.0);
            T vec_elem_j(0.0);
            atom_i->PolMatInteraction(*atom_j, mat_elem, vec_elem_i,
                                      vec_elem_j);
            
            pol_mat(i,j) = mat_elem;
            pol_mat(j,i) = mat_elem;
            vec_mat(i) += vec_elem_i;
            vec_mat(j) += vec_elem_j;
        }
    }
    
    for (unsigned i = 0; i < n_atoms; i++) {
        const Atom<T> *atom_i = atoms_molecules[i];
        
        T mat_elem(0.0);
        T vec_elem(0.0);
        atom_i->PolMatSelf(mat_elem, vec_elem);
        
        pol_mat(i,i) = mat_elem;
        vec_mat(i) += vec_elem;
        
        pol_mat(i,n_atoms) = atom_i->EffAtomicNumber();
        pol_mat(n_atoms,i) = atom_i->EffAtomicNumber();
        vec_mat(n_atoms) += atom_i->EffAtomicNumber();
    }
    
    Eigen::Matrix<T,Eigen::Dynamic,1> pol_coeffs;
    pol_coeffs = pol_mat.ldlt().solve(vec_mat);
    // pol_coeffs = pol_mat.lu().solve(vec_mat);
    
//    for (unsigned i = 0; i < n_atoms + 1; i++) {
//        for (unsigned j = 0; j < n_atoms + 1; j++) {
//            if (pol_mat(i,j) >= 0) {
//                char buffer[256];
//                snprintf(buffer,256," %8.3E  ",double(pol_mat(i,j)));
//                std::cout << buffer;
//            }
//            else {
//                char buffer[256];
//                snprintf(buffer,256,"%8.3E  ",double(pol_mat(i,j)));
//                std::cout << buffer;
//            }
//        }
//        char buffer[256];
//        if (vec_mat(i) >= 0) {
//            snprintf(buffer, 256, " %8.3E  ", double(vec_mat(i)));
//        }
//        else {
//            snprintf(buffer, 256, "%8.3E  ", double(vec_mat(i)));
//        }
//
//        std::cout << "|  " << buffer;
//        
//        if (pol_coeffs(i) >= 0) {
//            snprintf(buffer, 256, " %8.3E  ", double(pol_coeffs(i)));
//        }
//        else {
//            snprintf(buffer, 256, "%8.3E  ", double(pol_coeffs(i)));
//        }
//        
//        std::cout << buffer << std::endl;
//    }
//    
//    std::cout << std::endl;
    
    double aux_sum = 0;
    for (unsigned i = 0; i < n_atoms; i++) {
        atoms_molecules[i]->setPolCoeff(pol_coeffs(i));
    }
}

template<typename T>
T System<T>::SystemEnergy() {
    PolarizeMolecules();
    
    T energy = T(0.0);
    for (auto it1 = molecules.cbegin(); it1 != molecules.cend(); it1++) {
        for (auto it2 = it1 + 1; it2 != molecules.cend(); it2++) {
            energy += it1->InteractionEnergy(*it2);
        }
    }
    
    for (auto it1 = molecules.cbegin(); it1 != molecules.cend(); it1++) {
        energy += it1->SelfEnergy();
    }
    
    return energy;
}

template<typename T>
T System<T>::SystemInteractionEnergy() {
    T energy = SystemEnergy();
    
    for (std::string molec_name : names_molecules) {
        energy -= monomer_energies[molec_name];
    }
    
    return energy;
}

template<typename T>
void System<T>::MonomerPolarizeMolecules() {
    for (auto it = molecule_list.begin(); it != molecule_list.end(); it++) {
        it->second.Polarize();
    }
}

template<typename T>
void System<T>::getMonomerEnergies() {
    MonomerPolarizeMolecules();
    monomer_energies.clear();
    
    for (auto it = molecule_list.begin(); it != molecule_list.end(); it++) {
        monomer_energies[it->first] = double(it->second.SelfEnergy());
    }
}

template <typename T>
void System<T>::setECP() {
    for (auto it = molecules.begin(); it != molecules.end(); it++) {
        it->setECP();
    }
    
    for (auto it = molecule_list.begin(); it != molecule_list.end(); it++) {
        it->second.setECP();
    }
}

template <typename T>
const std::vector<std::string> System<T>::ListMoleculeTypes() const {
    std::vector<std::string> molec_names;
    molec_names.reserve(molecule_list.size());
    
    for (auto it = molecule_list.begin(); it != molecule_list.end(); it++) {
        molec_names.push_back(it->first);
    }
    
    return molec_names;
}

template <typename T>
const unsigned System<T>::MoleculeInstances(const std::string &molec_name) const {
    if (molecule_instances.find(molec_name) == molecule_instances.end()) {
        return 0;
    }
    
    return molecule_instances.find(molec_name)->second;
}

template <typename T>
const double System<T>::MoleculeMonomerEnergy(const std::string &molec_name) const {
    if (monomer_energies.find(molec_name) == monomer_energies.end()) {
        return 0;
    }
    
    return monomer_energies.find(molec_name)->second;
}

template <typename T>
void System<T>::setFullE() {
    for (auto it = molecules.begin(); it != molecules.end(); it++) {
        it->setFullE();
    }
    
    for (auto it = molecule_list.begin(); it != molecule_list.end(); it++) {
        it->second.setFullE();
    }
}

template <typename T>
void System<T>::setXCRule() {
    for (auto it = molecules.begin(); it != molecules.end(); it++) {
        it->setXCRule(xc_rule);
    }
    
    for (auto it = molecule_list.begin(); it != molecule_list.end(); it++) {
        it->second.setXCRule(xc_rule);
    }
}

template <typename T>
std::ostream& operator<<(std::ostream &os, const MolFFSim::System<T> &system) {
    os << "System Information" << std::endl;
    char buffer[MAX_PRINT_BUFFER_SIZE];
    for (int i = 0; i < 112; i++) {
        os << "-";
    }
    os << std::endl;
    
    snprintf(buffer, MAX_PRINT_BUFFER_SIZE,"%8s %25s %25s %25s %25s",
             "Element","x [Angstrom]","y [Angstrom]","z [Angstrom]",
             "Partial charge");
    os << buffer << std::endl;
    
    for (int i = 0; i < 112; i++) {
        os << "-";
    }
    os << std::endl;
    
    for (unsigned i = 0; i < system.Molecules().size(); i++) {
        if (!system.MoleculesToPrint()[i]) {
            continue;
        }
        
        os << system.Molecules()[i];
        
        for (int i = 0; i < 112; i++) {
            os << "-";
        }
        os << std::endl;
    }
    
    os << std::endl;
    os << std::endl;
    
    os << "Molecule/Fragment Information" << std::endl;
    
    for (int i = 0; i < 87; i++) {
        os << "-";
    }
    os << std::endl;
    
    snprintf(buffer, MAX_PRINT_BUFFER_SIZE,"%25s %25s %35s", "Name/Label",
             "Instances", "Total Energy [kJ/mol]");
    os << buffer << std::endl;
    
    for (int i = 0; i < 87; i++) {
        os << "-";
    }
    os << std::endl;
    
    for (std::string molec_name : system.ListMoleculeTypes()) {
        snprintf(buffer, MAX_PRINT_BUFFER_SIZE,"%25s %25u %35.5E", 
                 molec_name.c_str(), system.MoleculeInstances(molec_name),
                 system.MoleculeMonomerEnergy(molec_name));
        os << buffer << std::endl;
    }
    
    for (int i = 0; i < 87; i++) {
        os << "-";
    }
    
    os << std::endl;
    os << std::endl;
    
    return os;
}

// Explicit instantiation of all the types for the class System.
template class MolFFSim::System<double>;
template class MolFFSim::System<autodiff::var>;
template class MolFFSim::System<autodiff::dual>;

template std::ostream& operator<<(std::ostream &os, const MolFFSim::System<double> &system);
template std::ostream& operator<<(std::ostream &os, const MolFFSim::System<autodiff::var> &system);
template std::ostream& operator<<(std::ostream &os, const MolFFSim::System<autodiff::dual> &system);
