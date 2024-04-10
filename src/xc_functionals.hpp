#ifndef XC_FUNCTIONALS
#define XC_FUNCTIONALS

#include <cmath>

#include <autodiff/reverse/var.hpp>
#include <autodiff/forward/dual.hpp>

namespace MolFFSim {

// XC functional for the naive energy contribution.
template<typename T>
T NaiveModelE(const double lambda, const T &dist);

// XC functionals with spherically symmetric densities.
template<typename T>
T XCSpheSymm(const double lambda, const T &dist);

// XC functionals with a cylindrically symmetric density.
template<typename T>
T XCCylinSymm(const double lambda, const T &dist);

}

#endif
