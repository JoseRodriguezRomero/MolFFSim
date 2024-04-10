#define DEFAULT_ATOM_X              T(0.0)
#define DEFAULT_ATOM_Y              T(0.0)
#define DEFAULT_ATOM_Z              T(0.0)
#define DEFAULT_POL_COEFF           T(1.0)
#define DEFAULT_ATOMIC_NUMBER       1

#define DEFAULT_XC_A_COEFF          0.5
#define DEFAULT_XC_B_COEFF          0.5
#define DEFAULT_XC_C_COEFF          0.5
#define DEFAULT_XC_D_COEFF          0.5
#define DEFAULT_XC_A_POL_COEFF      0.5
#define DEFAULT_XC_B_POL_COEFF      0.5
#define DEFAULT_XC_C_POL_COEFF      0.5
#define DEFAULT_XC_D_POL_COEFF      0.5

#define DEFAULT_IS_PERIODIC         {false, false, false}
#define DEFAULT_PERIODIC_SIZES      {50, 50, 50}

#define DEFAULT_MOLECULE_POS_X      0.0
#define DEFAULT_MOLECULE_POS_Y      0.0
#define DEFAULT_MOLECULE_POS_Z      0.0
#define DEFAULT_MOLECULE_VEL_X      0.0
#define DEFAULT_MOLECULE_VEL_Y      0.0
#define DEFAULT_MOLECULE_VEL_Z      0.0
#define DEFAULT_MOLECULE_Q1         0.0
#define DEFAULT_MOLECULE_Q2         0.0
#define DEFAULT_MOLECULE_Q3         0.0
#define DEFAULT_MOLECULE_Q4         0.0
#define DEFAULT_MOLECULE_Q1_VEL     0.0
#define DEFAULT_MOLECULE_Q2_VEL     0.0
#define DEFAULT_MOLECULE_Q3_VEL     0.0
#define DEFAULT_MOLECULE_Q4_VEL     0.0

#define MAX_PRINT_BUFFER_SIZE       256

#define BOHR_TO_ANGSTROM            0.529177249

#include "molecule.hpp"

std::string labelFromAtomicNumber(const unsigned atomic_number) {
    switch (atomic_number) {
        case 1:
            return "H";
            break;
        case 2:
            return "He";
            break;
        case 3:
            return "Li";
            break;
        case 4:
            return "Be";
            break;
        case 5:
            return "B";
            break;
        case 6:
            return "C";
            break;
        case 7:
            return "N";
            break;
        case 8:
            return "O";
            break;
        case 9:
            return "F";
            break;
        case 10:
            return "Ne";
            break;
        case 11:
            return "Na";
            break;
        case 12:
            return "Mg";
            break;
        case 13:
            return "Al";
            break;
        case 14:
            return "Si";
            break;
        case 15:
            return "P";
            break;
        case 16:
            return "S";
            break;
        case 17:
            return "Cl";
            break;
        case 18:
            return "Ar";
            break;
        default:
            break;
    }

    // Dummy label for unsupported elements.
    return "XX";
}

void XCCombRule1(std::vector<double> &xc_ab, const std::vector<double> &xc_a,
                 const std::vector<double> &xc_b) {
    for (unsigned i = 0; i < 4; i++) {
        xc_ab[i] = (xc_a[i] * xc_b[i]) / (xc_a[i] + xc_b[i]);
    }
}

void XCCombRule2(std::vector<double> &xc_ab, const std::vector<double> &xc_a,
                 const std::vector<double> &xc_b) {
    for (unsigned i = 0; i < 4; i++) {
        xc_ab[i] = xc_a[i] + xc_b[i];
    }
}

void XCCombRule3(std::vector<double> &xc_ab, const std::vector<double> &xc_a,
                 const std::vector<double> &xc_b) {
    for (unsigned i = 0; i < 4; i++) {
        xc_ab[i] = xc_a[i] * xc_b[i];
    }
}

using namespace MolFFSim;

template<typename T>
Atom<T>::Atom() {
    pos.x() = DEFAULT_ATOM_X;
    pos.y() = DEFAULT_ATOM_Y;
    pos.z() = DEFAULT_ATOM_Z;
    atomic_number = DEFAULT_ATOMIC_NUMBER;
    
    pol_coeff = DEFAULT_POL_COEFF;
    
    xc_coeffs.reserve(4);
    xc_coeffs[0] = DEFAULT_XC_A_COEFF;
    xc_coeffs[1] = DEFAULT_XC_B_COEFF;
    xc_coeffs[2] = DEFAULT_XC_C_COEFF;
    xc_coeffs[3] = DEFAULT_XC_D_COEFF;
    
    is_periodic = DEFAULT_IS_PERIODIC;
    box_side_len = DEFAULT_PERIODIC_SIZES;
}

template<typename T>
Atom<T>::Atom(const Atom &other) {
    pos = other.Pos();
    atomic_number = other.AtomicNumber();
    
    pol_coeff = other.PolCoeff();
    cloud_c_coeffs = other.CloudCCoeffs();
    cloud_lambda_coeffs = other.CloudLambdaCoeffs();
    
    xc_coeffs.reserve(4);
    xc_coeffs[0] = XCCoeffA();
    xc_coeffs[1] = XCCoeffB();
    xc_coeffs[2] = XCCoeffC();
    xc_coeffs[3] = XCCoeffD();
    
    is_periodic = other.isPeriodic();
    box_side_len = other.PeriodicBoxSizes();
}

template<typename T>
Atom<T>::~Atom() {
}

template<typename T>
void Atom<T>::removeCloud(const unsigned cloud_index) {
    cloud_c_coeffs.erase(cloud_c_coeffs.begin() + cloud_index);
    cloud_lambda_coeffs.erase(cloud_lambda_coeffs.begin() + cloud_index);
}

template<typename T>
void Atom<T>::removeClouds(const unsigned from, const unsigned to) {
    cloud_c_coeffs.erase(cloud_c_coeffs.begin() + from,
                         cloud_c_coeffs.begin() + to + 1);
    cloud_lambda_coeffs.erase(cloud_lambda_coeffs.begin() + from,
                              cloud_lambda_coeffs.begin() + to + 1);
}

template<typename T>
void Atom<T>::addCloud(const double c_coeff, const double lambda_coeff) {
    cloud_c_coeffs.push_back(c_coeff);
    cloud_lambda_coeffs.push_back(lambda_coeff);
}

template<typename T>
T Atom<T>::SelfEnergy() const {
    T energy = 0;
    
    std::vector<double> comb_xc_coeffs;
    comb_xc_coeffs.reserve(4);
    
    switch (xc_rule) {
        case rule1:
            XCCombRule1(comb_xc_coeffs,xc_coeffs,xc_coeffs);
            break;
        case rule2:
            XCCombRule2(comb_xc_coeffs,xc_coeffs,xc_coeffs);
            break;
        default:
            XCCombRule3(comb_xc_coeffs,xc_coeffs,xc_coeffs);
            break;
    }
    
    // Cloud - Nuclei
    for (unsigned i = 0; i < cloud_c_coeffs.size(); i++) {
        double c_coeff = atomic_number*cloud_c_coeffs[i];
        double lambda_coeff = cloud_lambda_coeffs[i];
        
        energy -= c_coeff*NaiveModelE<T>(lambda_coeff,0);
        energy += c_coeff*comb_xc_coeffs[1]*XCSpheSymm<T>(lambda_coeff,0);
    }
    
    // Cloud - Cloud
    for (unsigned i = 0; i < cloud_c_coeffs.size(); i++) {
        for (unsigned j = i + 1; j < cloud_c_coeffs.size(); j++) {
            double c_coeff = cloud_c_coeffs[i] * cloud_c_coeffs[j];
            double lambda_coeff =
                (cloud_lambda_coeffs[i] * cloud_lambda_coeffs[j]) /
                (cloud_lambda_coeffs[i] + cloud_lambda_coeffs[j]);
            
            energy += c_coeff*NaiveModelE<T>(lambda_coeff,0);
            energy += c_coeff*comb_xc_coeffs[0]*XCSpheSymm<T>(lambda_coeff,0);
        }
    }
    
    return energy;
}

template<typename T>
T Atom<T>::InteractionEnergy(const Atom &other) const {
    T energy = 0;
    
    std::vector<double> comb_xc_coeffs;
    comb_xc_coeffs.reserve(4);
    
    switch (xc_rule) {
        case rule1:
            XCCombRule1(comb_xc_coeffs,xc_coeffs,other.XCCoeffs());
            break;
        case rule2:
            XCCombRule2(comb_xc_coeffs,xc_coeffs,other.XCCoeffs());
            break;
        default:
            XCCombRule3(comb_xc_coeffs,xc_coeffs,other.XCCoeffs());
            break;
    }
    
    T atoms_dist = Distance(other);
    
    // Nuclei - Nuclei
    energy += atomic_number*other.AtomicNumber()/atoms_dist;
    
    // Cloud - Nuclei
    for (unsigned i = 0; i < cloud_c_coeffs.size(); i++) {
        double c_coeff = other.AtomicNumber() * cloud_c_coeffs[i];
        double lambda_coeff = cloud_lambda_coeffs[i];
        
        energy -= c_coeff*NaiveModelE<T>(lambda_coeff,atoms_dist);
        energy +=
            c_coeff*comb_xc_coeffs[1]*XCSpheSymm<T>(lambda_coeff,atoms_dist);
        energy +=
            c_coeff*comb_xc_coeffs[3]*XCCylinSymm<T>(lambda_coeff,atoms_dist);
    }
    
    // Nuclei - Cloud
    for (unsigned i = 0; i < other.CloudCCoeffs().size(); i++) {
        double c_coeff = atomic_number * other.CloudCCoeffs()[i];
        double lambda_coeff = other.CloudLambdaCoeffs()[i];
        
        energy -= c_coeff*NaiveModelE<T>(lambda_coeff,atoms_dist);
        energy +=
            c_coeff*comb_xc_coeffs[1]*XCSpheSymm<T>(lambda_coeff,atoms_dist);
        energy +=
            c_coeff*comb_xc_coeffs[3]*XCCylinSymm<T>(lambda_coeff,atoms_dist);
    }
    
    // Cloud - Cloud
    for (unsigned i = 0; i < cloud_c_coeffs.size(); i++) {
        for (unsigned j = 0; j < other.CloudCCoeffs().size(); j++) {
            double c_coeff = cloud_c_coeffs[i] * other.CloudCCoeffs()[j];
            double lambda_coeff =
                (cloud_lambda_coeffs[i] * other.CloudLambdaCoeffs()[j]) /
                (cloud_lambda_coeffs[i] + other.CloudLambdaCoeffs()[j]);
            
            energy += c_coeff*NaiveModelE<T>(lambda_coeff,atoms_dist);
            energy += c_coeff*comb_xc_coeffs[0]*
                XCSpheSymm<T>(lambda_coeff,atoms_dist);
            energy += c_coeff*comb_xc_coeffs[2]*
                XCCylinSymm<T>(lambda_coeff,atoms_dist);
        }
    }
    
    return energy;
}

template<typename T>
void Atom<T>::operator=(const Atom& other) {
    pos = other.Pos();
    atomic_number = other.AtomicNumber();
    
    pol_coeff = other.PolCoeff();
    cloud_c_coeffs = other.CloudCCoeffs();
    cloud_lambda_coeffs = other.CloudLambdaCoeffs();
    
    xc_coeffs[0] = other.XCCoeffs()[0];
    xc_coeffs[1] = other.XCCoeffs()[1];
    xc_coeffs[2] = other.XCCoeffs()[2];
    xc_coeffs[3] = other.XCCoeffs()[3];
}

template<typename T>
bool Atom<T>::operator==(const Atom& other) const {
    if (this->Pos() != other.Pos()) {
        return false;
    }
    
    if (this->PolCoeff() != other.PolCoeff()) {
        return false;
    }
    
    if (this->CloudCCoeffs() != other.CloudCCoeffs()) {
        return false;
    }
    
    if (this->CloudLambdaCoeffs() != other.CloudLambdaCoeffs()) {
        return false;
    }
    
    if (this->XCCoeffA() != other.XCCoeffA()) {
        return false;
    }
    
    if (this->XCCoeffB() != other.XCCoeffB()) {
        return false;
    }
    
    if (this->XCCoeffC() != other.XCCoeffC()) {
        return false;
    }
    
    if (this->XCCoeffD() != other.XCCoeffD()) {
        return false;
    }
    
    return true;
}

template<typename T>
bool Atom<T>::operator!=(const Atom& other) const {
    return !((*this) == other);
}

template<typename T>
Molecule<T>::Molecule() {
    center_r.x() = DEFAULT_MOLECULE_POS_X;
    center_r.y() = DEFAULT_MOLECULE_POS_Y;
    center_r.z() = DEFAULT_MOLECULE_POS_Z;
    center_v.x() = DEFAULT_MOLECULE_VEL_X;
    center_v.y() = DEFAULT_MOLECULE_VEL_Y;
    center_v.z() = DEFAULT_MOLECULE_VEL_Z;
    
    center_qr.w() = DEFAULT_MOLECULE_Q1;
    center_qr.x() = DEFAULT_MOLECULE_Q2;
    center_qr.y() = DEFAULT_MOLECULE_Q3;
    center_qr.z() = DEFAULT_MOLECULE_Q4;
    center_qv.w() = DEFAULT_MOLECULE_Q1_VEL;
    center_qv.x() = DEFAULT_MOLECULE_Q2_VEL;
    center_qv.y() = DEFAULT_MOLECULE_Q3_VEL;
    center_qv.z() = DEFAULT_MOLECULE_Q4_VEL;
    
    is_periodic = DEFAULT_IS_PERIODIC;
    box_side_len = DEFAULT_PERIODIC_SIZES;
}

template<typename T>
Molecule<T>::Molecule(const Molecule& other) {
    center_r = other.Pos();
    center_v = other.Vel();
    center_qr = other.RotRotQ();
    center_qv = other.RotRotQVel();
    
    atoms = other.Atoms();
    atoms_rot = other.AtomsRotAndTrans();
    
    is_periodic = other.isPeriodic();
    box_side_len = other.PeriodicBoxSizes();
}

template<typename T>
Molecule<T>::~Molecule() {
}

template<typename T>
void Molecule<T>::addAtom(const MolFFSim::Atom<double> &atom) {
    atoms.push_back(atom);
    
    Eigen::Vector3<T> new_pos;
    new_pos.x() = atom.Pos().x();
    new_pos.y() = atom.Pos().y();
    new_pos.z() = atom.Pos().z();
    new_pos = center_qr*new_pos + center_r;
    
    Atom<T> new_atom;
    new_atom.setPos(new_pos);
    
    new_atom.setPolCoeff(atom.PolCoeff());
    new_atom.setAtomicNumber(atom.AtomicNumber());
    
    atoms_rot.push_back(new_atom);
}

template<typename T>
void Molecule<T>::removeAtom(const unsigned atom_index) {
    atoms.erase(atoms.begin() + atom_index);
    atoms_rot.erase(atoms_rot.begin() + atom_index);
}

template<typename T>
void Molecule<T>::removeAtoms(const unsigned from, const unsigned to) {
    atoms.erase(atoms.begin() + from, atoms.begin() + to + 1);
    atoms_rot.erase(atoms_rot.begin() + from, atoms_rot.begin() + to + 1);
}

template<typename T>
void Molecule<T>::addAtoms(const std::vector<MolFFSim::Atom<double>> &atoms) {
    for (auto it = atoms.begin(); it != atoms.end(); it++) {
        addAtom(*it);
    }
}

template<typename T>
void Molecule<T>::setAtomPolarizationCoeffs(const std::vector<T> pol_coeffs) {
    for (unsigned i = 0; i < pol_coeffs.size(); i++) {
        atoms_rot[i].setPolCoeff(pol_coeffs[i]);
    }
};

template<typename T>
void Molecule<T>::setPeriodic(const std::vector<bool>& is_periodic) {
    this->is_periodic = is_periodic;
    setPos(center_r);
    
    for (auto it = atoms_rot.begin(); it != atoms_rot.end(); it++) {
        it->setPeriodic(is_periodic);
    }
}

template<typename T>
void Molecule<T>::setPeriodicBoxSizes(const std::vector<double> box_side_len) {
    this->box_side_len = box_side_len;
    setPos(center_r);
    
    for (auto it = atoms_rot.begin(); it != atoms_rot.end(); it++) {
        it->setPeriodicBoxSizes(box_side_len);
    }
}

template<typename T>
void Molecule<T>::applyRotationAndTranslation() {
    for (unsigned i = 0; i < atoms.size(); i++) {
        Eigen::Vector3<T> atom_new_pos;
        atom_new_pos.x() = atoms[i].Pos().x();
        atom_new_pos.y() = atoms[i].Pos().y();
        atom_new_pos.z() = atoms[i].Pos().z();
        
        atom_new_pos = center_qr * atom_new_pos + center_r;
        atoms_rot[i].setPos(atom_new_pos);
    }
}

template<typename T>
T Molecule<T>::SelfEnergy() const {
    T energy = T(0);
    for (auto it1 = atoms_rot.begin(); it1 != atoms_rot.end(); it1++) {
        for (auto it2 = it1 + 1; it2 != atoms_rot.end(); it2++) {
            energy += it1->InteractionEnergy(*it2);
        }
        
        energy += it1->SelfEnergy();
    }
    
    return energy;
}

template<typename T>
T Molecule<T>::InteractionEnergy(const Molecule &other) const {
    T energy = T(0);
    for (auto it1 = atoms_rot.begin(); it1 != atoms_rot.end(); it1++) {
        auto aux_it_end = other.AtomsRotAndTrans().end();
        auto aux_it_begin = other.AtomsRotAndTrans().begin();
        for (auto it2 = aux_it_begin; it2 != aux_it_end; it2++) {
            energy += it1->InteractionEnergy(*it2);
        }
    }
    
    return energy;
}

template<typename T>
void Molecule<T>::setAnglesXYZ(const T &th_x, const T &th_y, const T &th_z) {
    const Eigen::Vector3<T> unit_x = Eigen::Vector3<T>::UnitX();
    const Eigen::Vector3<T> unit_y = Eigen::Vector3<T>::UnitY();
    const Eigen::Vector3<T> unit_z = Eigen::Vector3<T>::UnitZ();
    
    const Eigen::Matrix3<T> th_xr = 
        Eigen::AngleAxis<T>(th_x, unit_x).toRotationMatrix();
    const Eigen::Matrix3<T> th_yr = 
        Eigen::AngleAxis<T>(th_y, unit_y).toRotationMatrix();
    const Eigen::Matrix3<T> th_zr = 
        Eigen::AngleAxis<T>(th_z, unit_z).toRotationMatrix();
    
    Eigen::Matrix3<T> rot_mat = th_zr * th_yr * th_xr;
    center_qr = Eigen::Quaternion<T>(rot_mat);
}

template<typename T>
void Molecule<T>::operator=(const Molecule& other) {
    atoms = other.Atoms();
    atoms_rot = other.AtomsRotAndTrans();
    
    center_r = other.Pos();
    center_v = other.Vel();
    center_qr = other.RotRotQ();
    center_qv = other.RotRotQVel();
    
    is_periodic = other.isPeriodic();
    box_side_len = other.PeriodicBoxSizes();
    
    xc_rule = other.XCRule();
}

template<typename T>
bool Molecule<T>::operator==(const Molecule& other) const {
    return true;
}

template<typename T>
bool Molecule<T>::operator!=(const Molecule& other) const {
    return (*this == other);
}

template <typename T>
std::ostream& operator<<(std::ostream &os, const MolFFSim::Molecule<T> &molec) {
    char buffer[MAX_PRINT_BUFFER_SIZE];
    molec.SelfEnergy(); // Run this to ensure the molecule rotation has been
                        // updated before doing any printing.

//    for (int i = 0; i < 116; i++) {
//        os << "-";
//    }
//    os << std::endl;
//    
//    snprintf(buffer, MAX_PRINT_BUFFER_SIZE,"%12s %25s %25s %25s %25s",
//             "Element","x [Angstrom]","y [Angstrom]","z [Angstrom]",
//             "Partial charge");
//    os << buffer << std::endl;
//    
//    for (int i = 0; i < 116; i++) {
//        os << "-";
//    }
//    os << std::endl;
    
    // Print the data.
    auto it_end = molec.AtomsRotAndTrans().end();
    auto it_begin = molec.AtomsRotAndTrans().begin();
    for (auto it = it_begin; it != it_end; it++) {
        snprintf(buffer, MAX_PRINT_BUFFER_SIZE,
                 "%12s %25.8lf %25.8lf %25.8lf %25.8lf",
                 labelFromAtomicNumber(it->AtomicNumber()).c_str(),
                 double(it->Pos().x())*BOHR_TO_ANGSTROM,
                 double(it->Pos().y())*BOHR_TO_ANGSTROM,
                 double(it->Pos().z())*BOHR_TO_ANGSTROM,
                 double(it->ModelPartialCharge()));
        os << buffer << std::endl;
    }
    
//    for (int i = 0; i < 116; i++) {
//        os << "-";
//    }
//    os << std::endl;
    return os;
}

// Explicit instantiation of all the types for the class Atom.
template class MolFFSim::Atom<double>;
template class MolFFSim::Atom<autodiff::var>;
template class MolFFSim::Atom<autodiff::dual>;

// Explicit instantiation of all the types for the class Molecule.
template class MolFFSim::Molecule<double>;
template class MolFFSim::Molecule<autodiff::var>;
template class MolFFSim::Molecule<autodiff::dual>;

template std::ostream& operator<<(std::ostream &os, const MolFFSim::Molecule<double> &molec);
template std::ostream& operator<<(std::ostream &os, const MolFFSim::Molecule<autodiff::var> &molec);
template std::ostream& operator<<(std::ostream &os, const MolFFSim::Molecule<autodiff::dual> &molec);
