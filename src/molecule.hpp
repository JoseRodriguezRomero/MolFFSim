#ifndef CLASS_MOLECULE
#define CLASS_MOLECULE

#include <stdio.h>

#include <cmath>
#include <vector>
#include <string>
#include <algorithm>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <autodiff/reverse/var.hpp>
#include <autodiff/forward/dual.hpp>

#include "xc_functionals.hpp"

namespace MolFFSim {

template <typename T>
class Atom {
private:
    Eigen::Vector3<T> pos;                      // Atom data.
    unsigned atomic_number;                     // Atom data.
    
    T pol_coeff;                                // Electron cloud data.
    std::vector<double> cloud_c_coeffs;         // Electron cloud data.
    std::vector<double> cloud_lambda_coeffs;    // Electron cloud data.
    
    std::vector<double> xc_coeffs;
    
    enum XCRules {rule1, rule2, rule3};
    XCRules xc_rule = rule1;
    
    std::vector<bool> is_periodic;
    std::vector<double> box_side_len;
    
public:
    Atom();
    Atom(const Atom &other);
    ~Atom();
    
    inline Eigen::Vector3<T> Pos() const { return pos; }
    inline T PolCoeff() const { return pol_coeff; }
    inline unsigned AtomicNumber() const { return atomic_number; }
    
    inline T ModelPartialCharge() const {
        return atomic_number*(1.0 - pol_coeff);
    }
    
    inline const std::vector<bool>& isPeriodic() const { return is_periodic; }
    inline void setPeriodic(const std::vector<bool>& is_periodic) {
        this->is_periodic = is_periodic;
        setPos(pos);
    }
    
    const std::vector<double>& PeriodicBoxSizes() const {
        return box_side_len;
    }
    
    inline void setPeriodicBoxSizes(const std::vector<double> box_side_len) {
        this->box_side_len = box_side_len;
        setPos(pos);
    }
    
    inline T Distance(const Atom &other) const {
        Eigen::Vector3<T> dist_vec = pos - other.Pos();
        
        for (unsigned i = 0; i < 3; i++) {
            if (!is_periodic[i]) {
                continue;
            }
            
            while (dist_vec[i] > box_side_len[i] / 2.0) {
                dist_vec[i] -= box_side_len[i];
            }
            
            while (dist_vec[i] <= -box_side_len[i] / 2.0) {
                dist_vec[i] += box_side_len[i];
            }
        }
        
        return dist_vec.norm();
    }
    
    inline const std::vector<double>& CloudCCoeffs() const {
        return cloud_c_coeffs;
    }
    
    inline const std::vector<double>& CloudLambdaCoeffs() const {
        return cloud_lambda_coeffs;
    }
    
    inline const std::vector<double>& XCCoeffs() const { return xc_coeffs; }
    inline double XCCoeffA() const { return xc_coeffs[0]; }
    inline double XCCoeffB() const { return xc_coeffs[1]; }
    inline double XCCoeffC() const { return xc_coeffs[2]; }
    inline double XCCoeffD() const { return xc_coeffs[3]; }
    
    inline void setPos(const Eigen::Vector3<T> &pos) {
        this->pos = pos;
        
        for (unsigned i = 0; i < 3; i++) {
            if (!is_periodic[i]) {
                continue;
            }
            
            while (this->pos[i] > box_side_len[i] / 2.0) {
                this->pos[i] -= box_side_len[i];
            }
            
            while (this->pos[i] <= -box_side_len[i] / 2.0) {
                this->pos[i] += box_side_len[i];
            }
        }
    }
    
    inline void setPolCoeff(const T &pol_coeff) { this->pol_coeff = pol_coeff; }
    inline void setAtomicNumber(const unsigned atomic_number) {
        this->atomic_number = atomic_number;
    };
    
    void removeCloud(const unsigned cloud_index);
    void removeClouds(const unsigned from, const unsigned to);
    void addCloud(const double c_coeff, const double lambda_coeff);
    void addClouds(const std::vector<double> &c_coeffs,
                   const std::vector<double> &lambda_coeffs);
    
    T SelfEnergy() const;
    T InteractionEnergy(const Atom &other) const;
    
    void operator=(const Atom& other);
    bool operator==(const Atom& other) const;
    bool operator!=(const Atom& other) const;
};

template <typename T>
class Molecule {
private:
    std::vector<Atom<T>> atoms_rot;     // Totated and translated.
    std::vector<Atom<double>> atoms;    // Unrotated and untranslated.
    
    Eigen::Vector3<T> center_r;         // Position (center).
    Eigen::Vector3<T> center_v;         // Velocity (center).
    
    Eigen::Quaternion<T> center_qr;     // Quaternion rotation representation.
    Eigen::Quaternion<T> center_qv;     // Rate at which the quaternion
                                        // components change in time.
    
    enum XCRules {rule1, rule2, rule3};
    XCRules xc_rule = rule1;
    
    std::vector<bool> is_periodic;
    std::vector<double> box_side_len;
    
public:
    Molecule();
    Molecule(const Molecule &other);
    ~Molecule();
    
    inline Eigen::Vector3<T> Pos() const { return center_r; }
    inline Eigen::Vector3<T> Vel() const { return center_v; }
    inline Eigen::Quaternion<T> RotRotQ() const { return center_qr; }
    inline Eigen::Quaternion<T> RotRotQVel() const { return center_qv; }
    
    inline const std::vector<Atom<double>>& Atoms() const {
        return atoms;
    }
    
    inline const std::vector<Atom<T>>& AtomsRotAndTrans() const {
        return atoms_rot;
    }
    
    inline Molecule::XCRules XCRule() const { return xc_rule; }
    inline void setXCRule(Molecule::XCRules xc_rule) {
        this->xc_rule = xc_rule;
    }
    
    void removeAtom(const unsigned atom_index);
    void addAtom(const MolFFSim::Atom<double> &atom);
    void removeAtoms(const unsigned from, const unsigned to);
    void addAtoms(const std::vector<MolFFSim::Atom<double>> &atoms);
    
    void setAtomPolarizationCoeffs(const std::vector<T> pol_coeffs);
    
    void setPeriodic(const std::vector<bool>& is_periodic);
    inline const std::vector<bool>& isPeriodic() const { return is_periodic; }
    
    void setPeriodicBoxSizes(const std::vector<double> box_side_len);
    inline const std::vector<double>& PeriodicBoxSizes() const {
        return box_side_len;
    }
    
    inline void setPos(const Eigen::Vector3<T> &pos) {
        center_r = pos;
        
        for (unsigned i = 0; i < 3; i++) {
            if (!is_periodic[i]) {
                continue;
            }
            
            while (center_r[i] > box_side_len[i] / 2.0) {
                center_r[i] -= box_side_len[i];
            }
            
            while (center_r[i] <= -box_side_len[i] / 2.0) {
                center_r[i] += box_side_len[i];
            }
        }
    }
    
    inline void setVel(const Eigen::Vector3<T> &vel) { center_v = vel; }
    inline void setRotQ(const Eigen::Quaternion<T> &qr) { center_qr = qr; }
    inline void setRotQVel(const Eigen::Quaternion<T> &qv) { center_qv = qv; }
    void applyRotationAndTranslation(); // Invoke this function once the
                                        // rotations and translations of the
                                        // molecule are set.
    
    T SelfEnergy() const;
    T InteractionEnergy(const Molecule &other) const;
    
    void setAnglesXYZ(const T &th_x, const T &th_y, const T &th_z);
    
    void operator=(const Molecule& other);
    bool operator==(const Molecule& other) const;
    bool operator!=(const Molecule& other) const;
};

}

template <typename T>
std::ostream& operator<<(std::ostream &os, const MolFFSim::Molecule<T> &molec);

#endif
