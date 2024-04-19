#define DEFAULT_ATOM_X              T(0.0)
#define DEFAULT_ATOM_Y              T(0.0)
#define DEFAULT_ATOM_Z              T(0.0)
#define DEFAULT_POL_COEFF           T(1.0)
#define DEFAULT_ATOMIC_NUMBER       1

#define DEFAULT_XC_A_COEFF          0.0
#define DEFAULT_XC_B_COEFF          0.0
#define DEFAULT_XC_C_COEFF          0.0
#define DEFAULT_XC_D_COEFF          0.0
#define DEFAULT_XC_A_POL_COEFF      0.0
#define DEFAULT_XC_B_POL_COEFF      0.0
#define DEFAULT_XC_C_POL_COEFF      0.0
#define DEFAULT_XC_D_POL_COEFF      0.0

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

void XCCombRule1(std::vector<double> &xc_ab, const std::vector<double> &xc_a,
                 const std::vector<double> &xc_b) {
    for (unsigned i = 0; i < 8; i++) {
        xc_ab[i] = (xc_a[i] + xc_b[i]) / 2.0;
    }
}

void XCCombRule2(std::vector<double> &xc_ab, const std::vector<double> &xc_a,
                 const std::vector<double> &xc_b) {
    for (unsigned i = 0; i < 8; i++) {
        xc_ab[i] = sqrt(xc_a[i] * xc_b[i]);
    }
}

void XCCombRule3(std::vector<double> &xc_ab, const std::vector<double> &xc_a,
                 const std::vector<double> &xc_b) {
    for (unsigned i = 0; i < 8; i++) {
        xc_ab[i] = (xc_a[i] * xc_b[i]) / (xc_a[i] + xc_b[i]);
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
    
    xc_coeffs.reserve(8);
    xc_coeffs.push_back(DEFAULT_XC_A_COEFF);
    xc_coeffs.push_back(DEFAULT_XC_B_COEFF);
    xc_coeffs.push_back(DEFAULT_XC_C_COEFF);
    xc_coeffs.push_back(DEFAULT_XC_D_COEFF);
    xc_coeffs.push_back(DEFAULT_XC_A_POL_COEFF);
    xc_coeffs.push_back(DEFAULT_XC_B_POL_COEFF);
    xc_coeffs.push_back(DEFAULT_XC_C_POL_COEFF);
    xc_coeffs.push_back(DEFAULT_XC_D_POL_COEFF);
    
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
    xc_coeffs = other.XCCoeffs();
    
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
void Atom<T>::addClouds(const std::vector<double> &c_coeffs,
                        const std::vector<double> &lambda_coeffs) {
    cloud_c_coeffs.reserve(cloud_c_coeffs.size() + c_coeffs.size());
    cloud_lambda_coeffs.reserve(cloud_lambda_coeffs.size() + c_coeffs.size());
    
    for (unsigned i = 0; i < c_coeffs.size(); i++) {
        cloud_c_coeffs.push_back(c_coeffs[i]);
        cloud_lambda_coeffs.push_back(lambda_coeffs[i]);
    }
}

template<typename T>
T Atom<T>::SelfEnergy() const {
    T energy = 0;
    std::vector<double> comb_xc_coeffs = xc_coeffs;
    
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
        T c_coeff = z_eff*cloud_c_coeffs[i];
        double lambda_coeff = cloud_lambda_coeffs[i];
        
        if (i >= cloud_c_coeffs.size() - 3) {
            c_coeff *= pol_coeff;
        }
        
        energy -= c_coeff*NaiveModelE<T>(lambda_coeff,0);
        energy += c_coeff*comb_xc_coeffs[4]*XCSpheSymm<T>(lambda_coeff,0);
    }
    
    // Cloud - Cloud
    for (unsigned i = 0; i < cloud_c_coeffs.size(); i++) {
        for (unsigned j = 0; j < cloud_c_coeffs.size(); j++) {
            T c_coeff = cloud_c_coeffs[i] * cloud_c_coeffs[j];
            
            if (i >= cloud_c_coeffs.size() - 3) {
                c_coeff *= pol_coeff;
            }
            
            if (j >= cloud_c_coeffs.size() - 3) {
                c_coeff *= pol_coeff;
            }
            
            double lambda_coeff =
                (cloud_lambda_coeffs[i] * cloud_lambda_coeffs[j]) /
                (cloud_lambda_coeffs[i] + cloud_lambda_coeffs[j]);
            
            energy += c_coeff*NaiveModelE<T>(lambda_coeff,0);
            energy += c_coeff*comb_xc_coeffs[5]*XCSpheSymm<T>(lambda_coeff,0);
        }
    }
    
    return energy;
}

template<typename T>
T Atom<T>::InteractionEnergy(const Atom &other) const {
    T energy = 0;
    std::vector<double> comb_xc_coeffs = xc_coeffs;
    
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
    energy += z_eff*other.EffAtomicNumber()/atoms_dist;
    
    // Cloud - Nuclei
    for (unsigned i = 0; i < cloud_c_coeffs.size(); i++) {
        T c_coeff = other.EffAtomicNumber() * cloud_c_coeffs[i];
        double lambda_coeff = cloud_lambda_coeffs[i];
        
        if (i >= cloud_c_coeffs.size() - 3) {
            c_coeff *= pol_coeff;
        }
        
        energy -= c_coeff*NaiveModelE<T>(lambda_coeff,atoms_dist);
        energy +=
            c_coeff*comb_xc_coeffs[4]*XCSpheSymm<T>(lambda_coeff,atoms_dist);
        energy +=
            c_coeff*comb_xc_coeffs[6]*XCCylinSymm<T>(lambda_coeff,atoms_dist);
    }
    
    // Nuclei - Cloud
    for (unsigned i = 0; i < other.CloudCCoeffs().size(); i++) {
        T c_coeff = z_eff * other.CloudCCoeffs()[i];
        double lambda_coeff = other.CloudLambdaCoeffs()[i];
        
        if (i >= other.CloudCCoeffs().size() - 3) {
            c_coeff *= other.PolCoeff();
        }
        
        energy -= c_coeff*NaiveModelE<T>(lambda_coeff,atoms_dist);
        energy +=
            c_coeff*comb_xc_coeffs[4]*XCSpheSymm<T>(lambda_coeff,atoms_dist);
        energy +=
            c_coeff*comb_xc_coeffs[6]*XCCylinSymm<T>(lambda_coeff,atoms_dist);
    }
    
    // Cloud - Cloud
    for (unsigned i = 0; i < cloud_c_coeffs.size(); i++) {
        for (unsigned j = 0; j < other.CloudCCoeffs().size(); j++) {
            T c_coeff = cloud_c_coeffs[i] * other.CloudCCoeffs()[j];
            double lambda_coeff =
                (cloud_lambda_coeffs[i] * other.CloudLambdaCoeffs()[j]) /
                (cloud_lambda_coeffs[i] + other.CloudLambdaCoeffs()[j]);
            
            if (i >= cloud_c_coeffs.size() - 3) {
                c_coeff *= pol_coeff;
            }
            
            if (j >= other.CloudCCoeffs().size() - 3) {
                c_coeff *= other.PolCoeff();
            }
            
            energy += c_coeff*NaiveModelE<T>(lambda_coeff,atoms_dist);
            energy += c_coeff*comb_xc_coeffs[5]*
                XCSpheSymm<T>(lambda_coeff,atoms_dist);
            energy += c_coeff*comb_xc_coeffs[7]*
                XCCylinSymm<T>(lambda_coeff,atoms_dist);
        }
    }
    
    return energy;
}

template<typename T>
void Atom<T>::PolMatSelf(T &mat_elem, T &vec_elem) const {
    std::vector<double> comb_xc_coeffs = xc_coeffs;
    
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
        T c_coeff = z_eff*cloud_c_coeffs[i];
        double lambda_coeff = cloud_lambda_coeffs[i];
        
        if (i >= cloud_c_coeffs.size() - 3) {
            vec_elem += c_coeff*NaiveModelE<T>(lambda_coeff,0);
            vec_elem -= c_coeff*comb_xc_coeffs[0]*
                XCSpheSymm<T>(lambda_coeff,0);
        }
    }
        
    // Cloud - Cloud
    for (unsigned i = 0; i < cloud_c_coeffs.size(); i++) {
        for (unsigned j = 0; j < cloud_c_coeffs.size(); j++) {
            T c_coeff = cloud_c_coeffs[i] * cloud_c_coeffs[j];
            
            double lambda_coeff =
                (cloud_lambda_coeffs[i] * cloud_lambda_coeffs[j]) /
                (cloud_lambda_coeffs[i] + cloud_lambda_coeffs[j]);
            
            bool cond1 = (i >= (cloud_c_coeffs.size() - 3));
            bool cond2 = (j >= (cloud_c_coeffs.size() - 3));
            
            if (cond1 && cond2) {
                mat_elem += 2.0*c_coeff*NaiveModelE<T>(lambda_coeff,0);
                mat_elem += 2.0*c_coeff*comb_xc_coeffs[1]*
                    XCSpheSymm<T>(lambda_coeff,0);
            }
            else if (cond1 || cond2) {
                vec_elem -= c_coeff*NaiveModelE<T>(lambda_coeff,0);
                vec_elem -= c_coeff*comb_xc_coeffs[1]*
                    XCSpheSymm<T>(lambda_coeff,0);
            }
        }
    }
}

template<typename T>
void Atom<T>::PolMatInteraction(const Atom &other, T &mat_elem, 
                                T &vec_elem_i, T &vec_elem_j) const {
    
    std::vector<double> comb_xc_coeffs;
    comb_xc_coeffs.reserve(8);
    
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
    
    // Cloud - Nuclei
    for (unsigned i = 0; i < cloud_c_coeffs.size(); i++) {
        T c_coeff = other.EffAtomicNumber() * cloud_c_coeffs[i];
        double lambda_coeff = cloud_lambda_coeffs[i];
        
        if (i >= cloud_c_coeffs.size() - 3) {
            vec_elem_i += c_coeff*NaiveModelE<T>(lambda_coeff,atoms_dist);
            vec_elem_i -= c_coeff*comb_xc_coeffs[0]*
                XCSpheSymm<T>(lambda_coeff,atoms_dist);
            vec_elem_i -= c_coeff*comb_xc_coeffs[2]*
                XCCylinSymm<T>(lambda_coeff,atoms_dist);
        }
    }
    
    // Nuclei - Cloud
    for (unsigned i = 0; i < other.CloudCCoeffs().size(); i++) {
        T c_coeff = z_eff * other.CloudCCoeffs()[i];
        double lambda_coeff = other.CloudLambdaCoeffs()[i];
        
        if (i >= other.CloudCCoeffs().size() - 3) {
            vec_elem_j += c_coeff*NaiveModelE<T>(lambda_coeff,atoms_dist);
            vec_elem_j -= c_coeff*comb_xc_coeffs[0]*
                XCSpheSymm<T>(lambda_coeff,atoms_dist);
            vec_elem_j -= c_coeff*comb_xc_coeffs[2]*
                XCCylinSymm<T>(lambda_coeff,atoms_dist);
        }
    }
    
    // Cloud - Cloud
    for (unsigned i = 0; i < cloud_c_coeffs.size(); i++) {
        for (unsigned j = 0; j < other.CloudCCoeffs().size(); j++) {
            T c_coeff = cloud_c_coeffs[i] * other.CloudCCoeffs()[j];
            double lambda_coeff =
                (cloud_lambda_coeffs[i] * other.CloudLambdaCoeffs()[j]) /
                (cloud_lambda_coeffs[i] + other.CloudLambdaCoeffs()[j]);
            
            bool cond1 = (i >= cloud_c_coeffs.size() - 3);
            bool cond2 = (j >= other.CloudCCoeffs().size() - 3);
            
            if (cond1 && cond2) {
                mat_elem += c_coeff*NaiveModelE<T>(lambda_coeff,atoms_dist);
                mat_elem += c_coeff*comb_xc_coeffs[1]*
                    XCSpheSymm<T>(lambda_coeff,atoms_dist);
                mat_elem += c_coeff*comb_xc_coeffs[3]*
                    XCCylinSymm<T>(lambda_coeff,atoms_dist);
            }
            else if (cond1) {
                vec_elem_i -= c_coeff*NaiveModelE<T>(lambda_coeff,atoms_dist);
                vec_elem_i -= c_coeff*comb_xc_coeffs[1]*
                    XCSpheSymm<T>(lambda_coeff,atoms_dist);
                vec_elem_i -= c_coeff*comb_xc_coeffs[3]*
                    XCCylinSymm<T>(lambda_coeff,atoms_dist);
            }
            else if (cond2) {
                vec_elem_j -= c_coeff*NaiveModelE<T>(lambda_coeff,atoms_dist);
                vec_elem_j -= c_coeff*comb_xc_coeffs[1]*
                    XCSpheSymm<T>(lambda_coeff,atoms_dist);
                vec_elem_j -= c_coeff*comb_xc_coeffs[3]*
                    XCCylinSymm<T>(lambda_coeff,atoms_dist);
            }
        }
    }
}

template<typename T>
void Atom<T>::operator=(const Atom& other) {
    pos = other.Pos();
    atomic_number = other.AtomicNumber();
    
    pol_coeff = other.PolCoeff();
    cloud_c_coeffs = other.CloudCCoeffs();
    cloud_lambda_coeffs = other.CloudLambdaCoeffs();
    xc_coeffs = other.XCCoeffs();
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
    center_qr = other.RotQ();
    center_qv = other.RotQVel();
    
    atoms = other.const_Atoms();
    atoms_rot = other.const_AtomsRotAndTrans();
    
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
void Molecule<T>::setXCRule(MolFFSim::XCRules xc_rule) {    
    for (auto it = atoms.begin(); it != atoms.end(); it++) {
        it->setXCRule(xc_rule);
    }
    
    for (auto it = atoms_rot.begin(); it != atoms_rot.end(); it++) {
        it->setXCRule(xc_rule);
    }
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
void Molecule<T>::createElectronClouds(
    const std::unordered_map<unsigned,std::vector<double>> &c_coeffs,
    const std::unordered_map<unsigned,std::vector<double>> &lambda_coeffs) {
    for (auto it = atoms_rot.begin(); it != atoms_rot.end(); it++) {
        it->addClouds(c_coeffs.find(it->AtomicNumber())->second,
                      lambda_coeffs.find(it->AtomicNumber())->second);
    }
    
    for (auto it = atoms.begin(); it != atoms.end(); it++) {
        it->addClouds(c_coeffs.find(it->AtomicNumber())->second,
                      lambda_coeffs.find(it->AtomicNumber())->second);
    }
}

template<typename T>
void Molecule<T>::setXCCoefficients(
    const std::unordered_map<unsigned,std::vector<double>> &xc_coeffs) {
    for (auto it = atoms_rot.begin(); it != atoms_rot.end(); it++) {
        it->setXCCoeffs(xc_coeffs.find(it->AtomicNumber())->second);
    }
    
    for (auto it = atoms.begin(); it != atoms.end(); it++) {
        it->setXCCoeffs(xc_coeffs.find(it->AtomicNumber())->second);
    }
}

template<typename T>
void Molecule<T>::setECP() {
    for (auto it = atoms_rot.begin(); it != atoms_rot.end(); it++) {
        it->setECP();
    }
    
    for (auto it = atoms.begin(); it != atoms.end(); it++) {
        it->setECP();
    }
}

template<typename T>
void Molecule<T>::setFullE() {
    for (auto it = atoms_rot.begin(); it != atoms_rot.end(); it++) {
        it->setFullE();
    }
    
    for (auto it = atoms.begin(); it != atoms.end(); it++) {
        it->setFullE();
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
        auto aux_it_end = other.const_AtomsRotAndTrans().cend();
        auto aux_it_begin = other.const_AtomsRotAndTrans().cbegin();
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
void Molecule<T>::Polarize() {
    unsigned n_atoms = atoms_rot.size();
    Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> pol_mat(n_atoms+1,n_atoms+1);
    Eigen::Matrix<T,Eigen::Dynamic,1> vec_mat(n_atoms+1);
    
    pol_mat.setZero();
    vec_mat.setZero();
    
    for (unsigned i = 0; i < n_atoms; i++) {
        T mat_elem(0.0);
        T vec_elem(0.0);
        atoms_rot[i].PolMatSelf(mat_elem, vec_elem);
        
        pol_mat(i,i) = mat_elem;
        vec_mat(i) += vec_elem;
        
        pol_mat(i,n_atoms) =
            ECPEffectiveAtomicNumber(atoms_rot[i].AtomicNumber());
        pol_mat(n_atoms,i) = 
            ECPEffectiveAtomicNumber(atoms_rot[i].AtomicNumber());
        vec_mat(n_atoms) += 
            ECPEffectiveAtomicNumber(atoms_rot[i].AtomicNumber());
    }
    
    for (unsigned i = 0; i < n_atoms; i++) {
        for (unsigned j = i + 1; j < n_atoms; j++) {
            T mat_elem(0.0);
            T vec_elem_i(0.0);
            T vec_elem_j(0.0);
            atoms_rot[i].PolMatInteraction(atoms_rot[j], mat_elem, vec_elem_i,
                                           vec_elem_j);
            
             vec_mat(i) += vec_elem_i;
             vec_mat(j) += vec_elem_j;
             pol_mat(i,j) = mat_elem;
             pol_mat(j,i) = mat_elem;
        }
    }
    
    Eigen::Matrix<T,Eigen::Dynamic,1> pol_coeffs;
    pol_coeffs = pol_mat.ldlt().solve(vec_mat);
    
    for (unsigned i = 0; i < n_atoms; i++) {
        atoms_rot[i].setPolCoeff(pol_coeffs(i));
    }
}

template<typename T>
void Molecule<T>::operator=(const Molecule& other) {
    atoms = other.const_Atoms();
    atoms_rot = other.const_AtomsRotAndTrans();
    
    center_r = other.Pos();
    center_v = other.Vel();
    center_qr = other.RotQ();
    center_qv = other.RotQVel();
    
    is_periodic = other.isPeriodic();
    box_side_len = other.PeriodicBoxSizes();
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

//    for (int i = 0; i < 112; i++) {
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
    auto it_end = molec.const_AtomsRotAndTrans().cend();
    auto it_begin = molec.const_AtomsRotAndTrans().cbegin();
    for (auto it = it_begin; it != it_end; it++) {
        snprintf(buffer, MAX_PRINT_BUFFER_SIZE,
                 "%8s %25.5lf %25.5lf %25.5lf %25.5lf",
                 LabelFromAtomicNumber(it->AtomicNumber()).c_str(),
                 double(it->Pos().x())*BOHR_TO_ANGSTROM,
                 double(it->Pos().y())*BOHR_TO_ANGSTROM,
                 double(it->Pos().z())*BOHR_TO_ANGSTROM,
                 double(it->ModelPartialCharge()));
        os << buffer << std::endl;
    }
    
//    for (int i = 0; i < 112; i++) {
//        os << "-";
//    }
//    os << std::endl;
    return os;
}

// Explicit instantiation of all the types for the class Atom.
template class MolFFSim::Atom<double>;
template class MolFFSim::Atom<autodiff::dual>;

// Explicit instantiation of all the types for the class Molecule.
template class MolFFSim::Molecule<double>;
template class MolFFSim::Molecule<autodiff::dual>;

template std::ostream& operator<<(std::ostream &os, const MolFFSim::Molecule<double> &molec);
template std::ostream& operator<<(std::ostream &os, const MolFFSim::Molecule<autodiff::dual> &molec);
