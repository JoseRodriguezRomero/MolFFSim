#define M_PI_SQRT                   1.7724538509055159
#define MIN_ATOM_POTENTIAL_DIST     1.0E-6

#define POW_TABLE_SIZE              21

#define FACTORIAL_2                 2.0
#define FACTORIAL_3                 6.0
#define FACTORIAL_4                 24.0
#define FACTORIAL_5                 120.0
#define FACTORIAL_6                 720.0
#define FACTORIAL_7                 5040.0
#define FACTORIAL_8                 40320.0
#define FACTORIAL_9                 362880.0
#define FACTORIAL_10                3628800.0
#define FACTORIAL_11                39916800.0
#define FACTORIAL_12                479001600.0
#define FACTORIAL_13                6227020800.0
#define FACTORIAL_14                87178291200.0
#define FACTORIAL_15                1307674368000.0
#define FACTORIAL_16                20922789888000.0
#define FACTORIAL_17                355687428096000.0
#define FACTORIAL_18                6402373705728000.0
#define FACTORIAL_19                121645100408832000.0
#define FACTORIAL_20                2432902008176640000.0
#define FACTORIAL_21                51090942171709440000.0

#define APPROX_RADIUS               25.0

#include "xc_functionals.hpp"

template<typename T>
T XCEnergy1(const std::vector<double> &lambda_pows, 
            const std::vector<T> &dist_pows) {
    // Multiply this with sqrt(lambda) * exp(-lambda*dist^2) / M_PI_SQRT
    return -(4*lambda_pows[1]);
}

template<typename T>
T XCEnergy2(const std::vector<double> &lambda_pows, 
            const std::vector<T> &dist_pows) {
    // Multiply this with sqrt(lambda) * exp(-lambda*dist^2) / M_PI_SQRT
    return -(16*dist_pows[2]*lambda_pows[3]-24*lambda_pows[2]);
}

template<typename T>
T XCEnergy3(const std::vector<double> &lambda_pows, 
            const std::vector<T> &dist_pows) {
    // Multiply this with sqrt(lambda) * exp(-lambda*dist^2) / M_PI_SQRT
    return -(64*dist_pows[4]*lambda_pows[5]-320*dist_pows[2]*lambda_pows[4]+
             240*lambda_pows[3]);
}

template<typename T>
T XCEnergy4(const std::vector<double> &lambda_pows, 
            const std::vector<T> &dist_pows) {
    // Multiply this with sqrt(lambda) * exp(-lambda*dist^2) / M_PI_SQRT
    return -(256*dist_pows[6]*lambda_pows[7]-2688*dist_pows[4]*lambda_pows[6]+
             6720*dist_pows[2]*lambda_pows[5]-3360*lambda_pows[4]);
}

template<typename T>
T XCEnergy5(const std::vector<double> &lambda_pows, 
            const std::vector<T> &dist_pows) {
    // Multiply this with sqrt(lambda) * exp(-lambda*dist^2) / M_PI_SQRT
    return -(1024*dist_pows[8]*lambda_pows[9]-18432*dist_pows[6]*lambda_pows[8]+
             96768*dist_pows[4]*lambda_pows[7]-161280*dist_pows[2]*
             lambda_pows[6]+60480*lambda_pows[5]);
}

template<typename T>
T XCEnergy6(const std::vector<double> &lambda_pows, 
            const std::vector<T> &dist_pows) {
    // Multiply this with sqrt(lambda) * exp(-lambda*dist^2) / M_PI_SQRT
    return -(4096*dist_pows[10]*lambda_pows[11]-112640*dist_pows[8]*
             lambda_pows[10]+1013760*dist_pows[6]*lambda_pows[9]-3548160*
             dist_pows[4]*lambda_pows[8]+4435200*dist_pows[2]*lambda_pows[7]-
             1330560*lambda_pows[6]);
}

template<typename T>
T XCEnergy7(const std::vector<double> &lambda_pows, 
            const std::vector<T> &dist_pows) {
    // Multiply this with sqrt(lambda) * exp(-lambda*dist^2) / M_PI_SQRT
    return -(16384*dist_pows[12]*lambda_pows[13]-638976*dist_pows[10]*
             lambda_pows[12]+8785920*dist_pows[8]*lambda_pows[11]-52715520*
             dist_pows[6]*lambda_pows[10]+138378240*dist_pows[4]*lambda_pows[9]-
             138378240*dist_pows[2]*lambda_pows[8]+34594560*lambda_pows[7]);
}

template<typename T>
T XCEnergy8(const std::vector<double> &lambda_pows, 
            const std::vector<T> &dist_pows) {
    // Multiply this with sqrt(lambda) * exp(-lambda*dist^2) / M_PI_SQRT
    return -(65536*dist_pows[14]*lambda_pows[15]-3440640*dist_pows[12]*
             lambda_pows[14]+67092480*dist_pows[10]*lambda_pows[13]-615014400*
             dist_pows[8]*lambda_pows[12]+2767564800*dist_pows[6]*
             lambda_pows[11]-5811886080*dist_pows[4]*lambda_pows[10]+4843238400*
             dist_pows[2]*lambda_pows[9]-1037836800*lambda_pows[8]);
}

template<typename T>
T XCEnergy9(const std::vector<double> &lambda_pows, 
            const std::vector<T> &dist_pows) {
    // Multiply this with sqrt(lambda) * exp(-lambda*dist^2) / M_PI_SQRT
    return -(262144*dist_pows[16]*lambda_pows[17]-17825792*dist_pows[14]*
             lambda_pows[16]+467927040*dist_pows[12]*lambda_pows[15]-6083051520*
             dist_pows[10]*lambda_pows[14]+41820979200*dist_pows[8]*
             lambda_pows[13]-150555525120*dist_pows[6]*lambda_pows[12]+
             263472168960*dist_pows[4]*lambda_pows[11]-188194406400*
             dist_pows[2]*lambda_pows[10]+35286451200*lambda_pows[9]);
}

template<typename T>
T XCEnergy10(const std::vector<double> &lambda_pows, 
             const std::vector<T> &dist_pows) {
    // Multiply this with sqrt(lambda) * exp(-lambda*dist^2) / M_PI_SQRT
    return -(1048576*dist_pows[18]*lambda_pows[19]-89653248*dist_pows[16]*
             lambda_pows[18]+3048210432*dist_pows[14]*lambda_pows[17]-
             53343682560*dist_pows[12]*lambda_pows[16]+520100904960*
             dist_pows[10]*lambda_pows[15]-2860554977280*dist_pows[8]*
             lambda_pows[14]+8581664931840*dist_pows[6]*lambda_pows[13]-
             12872497397760*dist_pows[4]*lambda_pows[12]+8045310873600*
             dist_pows[2]*lambda_pows[11]-1340885145600*lambda_pows[10]);
}

template<typename T>
T XCEnergyD1(const std::vector<double> &lambda_pows, 
             const std::vector<T> &dist_pows) {
    // Multiply this with sqrt(lambda) * exp(-lambda*dist^2) / M_PI_SQRT
    return 8*dist_pows[1]*lambda_pows[2];
}

template<typename T>
T XCEnergyD2(const std::vector<double> &lambda_pows, 
             const std::vector<T> &dist_pows) {
    // Multiply this with sqrt(lambda) * exp(-lambda*dist^2) / M_PI_SQRT
    return 32*dist_pows[3]*lambda_pows[4]-80*dist_pows[1]*lambda_pows[3];
}

template<typename T>
T XCEnergyD3(const std::vector<double> &lambda_pows, 
             const std::vector<T> &dist_pows) {
    // Multiply this with sqrt(lambda) * exp(-lambda*dist^2) / M_PI_SQRT
    return 128*dist_pows[5]*lambda_pows[6]-896*dist_pows[3]*lambda_pows[5]+
        1120*dist_pows[1]*lambda_pows[4];
}

template<typename T>
T XCEnergyD4(const std::vector<double> &lambda_pows, 
             const std::vector<T> &dist_pows) {
    // Multiply this with sqrt(lambda) * exp(-lambda*dist^2) / M_PI_SQRT
    return 512*dist_pows[7]*lambda_pows[8]-6912*dist_pows[5]*lambda_pows[7]+
        24192*dist_pows[3]*lambda_pows[6]-20160*dist_pows[1]*lambda_pows[5];
}

template<typename T>
T XCEnergyD5(const std::vector<double> &lambda_pows, 
             const std::vector<T> &dist_pows) {
    // Multiply this with sqrt(lambda) * exp(-lambda*dist^2) / M_PI_SQRT
    return 2048*dist_pows[9]*lambda_pows[10]-45056*dist_pows[7]*lambda_pows[9]+
        304128*dist_pows[5]*lambda_pows[8]-709632*dist_pows[3]*lambda_pows[7]+
        443520*dist_pows[1]*lambda_pows[6];
}

template<typename T>
T XCEnergyD6(const std::vector<double> &lambda_pows, 
             const std::vector<T> &dist_pows) {
    // Multiply this with sqrt(lambda) * exp(-lambda*dist^2) / M_PI_SQRT
    return 8192*dist_pows[11]*lambda_pows[12]-266240*dist_pows[9]*
        lambda_pows[11]+2928640*dist_pows[7]*lambda_pows[10]-13178880*
        dist_pows[5]*lambda_pows[9]+23063040*dist_pows[3]*lambda_pows[8]-
        11531520*dist_pows[1]*lambda_pows[7];
}

template<typename T>
T XCEnergyD7(const std::vector<double> &lambda_pows, 
             const std::vector<T> &dist_pows) {
    // Multiply this with sqrt(lambda) * exp(-lambda*dist^2) / M_PI_SQRT
    return 32768*dist_pows[13]*lambda_pows[14]-1474560*dist_pows[11]*
        lambda_pows[13]+23961600*dist_pows[9]*lambda_pows[12]-175718400*
        dist_pows[7]*lambda_pows[11]+593049600*dist_pows[5]*lambda_pows[10]-
        830269440*dist_pows[3]*lambda_pows[9]+345945600*dist_pows[1]*
        lambda_pows[8];
}

template<typename T>
T XCEnergyD8(const std::vector<double> &lambda_pows, 
             const std::vector<T> &dist_pows) {
    // Multiply this with sqrt(lambda) * exp(-lambda*dist^2) / M_PI_SQRT
    return 131072*dist_pows[15]*lambda_pows[16]-7798784*dist_pows[13]*
        lambda_pows[15]+175472640*dist_pows[11]*lambda_pows[14]-1900953600*
        dist_pows[9]*lambda_pows[13]+10455244800*dist_pows[7]*lambda_pows[12]-
        28229160960*dist_pows[5]*lambda_pows[11]+32934021120*dist_pows[3]*
        lambda_pows[10]-11762150400*dist_pows[1]*lambda_pows[9];
}

template<typename T>
T XCEnergyD9(const std::vector<double> &lambda_pows, 
             const std::vector<T> &dist_pows) {
    // Multiply this with sqrt(lambda) * exp(-lambda*dist^2) / M_PI_SQRT
    return 524288*dist_pows[17]*lambda_pows[18]-39845888*dist_pows[15]*
        lambda_pows[17]+1185415168*dist_pows[13]*lambda_pows[16]-17781227520*
        dist_pows[11]*lambda_pows[15]+144472473600*dist_pows[9]*lambda_pows[14]-
        635678883840*dist_pows[7]*lambda_pows[13]+1430277488640*dist_pows[5]*
        lambda_pows[12]-1430277488640*dist_pows[3]*lambda_pows[11]+446961715200*
        dist_pows[1]*lambda_pows[10];
}

template<typename T>
T XCEnergyD10(const std::vector<double> &lambda_pows, 
              const std::vector<T> &dist_pows) {
    // Multiply this with sqrt(lambda) * exp(-lambda*dist^2) / M_PI_SQRT
    return 2097152*dist_pows[19]*lambda_pows[20]-198180864*dist_pows[17]*
        lambda_pows[19]+7530872832*dist_pows[15]*lambda_pows[18]-149362311168*
        dist_pows[13]*lambda_pows[17]+1680326000640*dist_pows[11]*
        lambda_pows[16]-10922119004160*dist_pows[9]*lambda_pows[15]+
        40047769681920*dist_pows[7]*lambda_pows[14]-77234984386560*dist_pows[5]*
        lambda_pows[13]+67580611338240*dist_pows[3]*lambda_pows[12]-
        18772392038400*dist_pows[1]*lambda_pows[11];
}

template<typename T>
T MolFFSim::NaiveModelE(const double lambda, const T &dist) {
    if (dist <= MIN_ATOM_POTENTIAL_DIST) {
        return 2.0*sqrt(lambda) / M_PI_SQRT;
    }
        
    if ((lambda*dist*dist) >= APPROX_RADIUS) {
        return 1.0 / dist;
    }
    
    return erf(sqrt(lambda)*dist) / dist;
}

template<typename T>
T MolFFSim::XCSpheSymm(const double lambda, const T &dist) {
    if ((lambda*dist*dist) >= APPROX_RADIUS) {
        return T(0.0);
    }
    
    std::vector<T> dist_pows;
    std::vector<double> lambda_pows;
    
    dist_pows.reserve(POW_TABLE_SIZE);
    lambda_pows.reserve(POW_TABLE_SIZE);
    
    if (dist == 0) {
        for (unsigned i = 0; i < POW_TABLE_SIZE; i++) {
            dist_pows[i] = 0;
            lambda_pows[i] = pow(lambda,i);
        }
    }
    else {
        for (unsigned i = 0; i < POW_TABLE_SIZE; i++) {
            dist_pows[i] = pow(dist,i);
            lambda_pows[i] = pow(lambda,i);
        }
    }
    
    // THIS IS PROBLEMATIC WITH AUTODIFF!!!
    // T aux_c = sqrt(lambda) /  (exp(lambda*dist_pows[2]) * M_PI_SQRT);
    T aux_c = sqrt(lambda) * exp(-lambda*dist_pows[2]) / M_PI_SQRT;
    
    T xc_energy = T(0.0);
    xc_energy -= XCEnergy1(lambda_pows, dist_pows) / FACTORIAL_2;
    xc_energy += XCEnergy2(lambda_pows, dist_pows) / FACTORIAL_4;
    xc_energy -= XCEnergy3(lambda_pows, dist_pows) / FACTORIAL_6;
    xc_energy += XCEnergy4(lambda_pows, dist_pows) / FACTORIAL_8;
    xc_energy -= XCEnergy5(lambda_pows, dist_pows) / FACTORIAL_10;
    xc_energy += XCEnergy6(lambda_pows, dist_pows) / FACTORIAL_12;
    xc_energy -= XCEnergy7(lambda_pows, dist_pows) / FACTORIAL_14;
    xc_energy += XCEnergy8(lambda_pows, dist_pows) / FACTORIAL_16;
    xc_energy -= XCEnergy9(lambda_pows, dist_pows) / FACTORIAL_18;
    xc_energy += XCEnergy10(lambda_pows, dist_pows) / FACTORIAL_20;
    
    return aux_c * xc_energy;
}

template<typename T>
T MolFFSim::XCCylinSymm(const double lambda, const T &dist) {
    if ((lambda*dist*dist) >= APPROX_RADIUS) {
        return T(0.0);
    }
    
    std::vector<T> dist_pows;
    std::vector<double> lambda_pows;
    
    dist_pows.reserve(POW_TABLE_SIZE);
    lambda_pows.reserve(POW_TABLE_SIZE);
    
    if (dist == 0) {
        for (unsigned i = 0; i < POW_TABLE_SIZE; i++) {
            dist_pows[i] = 0;
            lambda_pows[i] = pow(lambda,i);
        }
    }
    else {
        for (unsigned i = 0; i < POW_TABLE_SIZE; i++) {
            dist_pows[i] = pow(dist,i);
            lambda_pows[i] = pow(lambda,i);
        }
    }
    
    // THIS IS PROBLEMATIC WITH AUTODIFF!!!
    // T aux_c = sqrt(lambda) /  (exp(lambda*dist_pows[2]) * M_PI_SQRT);
    T aux_c = sqrt(lambda) * exp(-lambda*dist_pows[2]) / M_PI_SQRT;
    
    T xc_energy = T(0.0);
    xc_energy += XCEnergyD1(lambda_pows, dist_pows) / FACTORIAL_3;
    xc_energy -= XCEnergyD2(lambda_pows, dist_pows) / FACTORIAL_5;
    xc_energy += XCEnergyD3(lambda_pows, dist_pows) / FACTORIAL_7;
    xc_energy -= XCEnergyD4(lambda_pows, dist_pows) / FACTORIAL_9;
    xc_energy += XCEnergyD5(lambda_pows, dist_pows) / FACTORIAL_11;
    xc_energy -= XCEnergyD6(lambda_pows, dist_pows) / FACTORIAL_13;
    xc_energy += XCEnergyD7(lambda_pows, dist_pows) / FACTORIAL_15;
    xc_energy -= XCEnergyD8(lambda_pows, dist_pows) / FACTORIAL_17;
    xc_energy += XCEnergyD9(lambda_pows, dist_pows) / FACTORIAL_19;
    xc_energy -= XCEnergyD10(lambda_pows, dist_pows) / FACTORIAL_21;
    
    return aux_c * xc_energy;
}

// Explicit instantiation of all the auxiliary types for the XC functionals.
template double XCEnergy1(const std::vector<double> &lambda_pows, const std::vector<double> &dist_pows);
template double XCEnergy2(const std::vector<double> &lambda_pows, const std::vector<double> &dist_pows);
template double XCEnergy3(const std::vector<double> &lambda_pows, const std::vector<double> &dist_pows);
template double XCEnergy4(const std::vector<double> &lambda_pows, const std::vector<double> &dist_pows);
template double XCEnergy5(const std::vector<double> &lambda_pows, const std::vector<double> &dist_pows);
template double XCEnergy6(const std::vector<double> &lambda_pows, const std::vector<double> &dist_pows);
template double XCEnergy7(const std::vector<double> &lambda_pows, const std::vector<double> &dist_pows);
template double XCEnergy8(const std::vector<double> &lambda_pows, const std::vector<double> &dist_pows);
template double XCEnergy9(const std::vector<double> &lambda_pows, const std::vector<double> &dist_pows);
template double XCEnergy10(const std::vector<double> &lambda_pows, const std::vector<double> &dist_pows);

template double XCEnergyD1(const std::vector<double> &lambda_pows, const std::vector<double> &dist_pows);
template double XCEnergyD2(const std::vector<double> &lambda_pows, const std::vector<double> &dist_pows);
template double XCEnergyD3(const std::vector<double> &lambda_pows, const std::vector<double> &dist_pows);
template double XCEnergyD4(const std::vector<double> &lambda_pows, const std::vector<double> &dist_pows);
template double XCEnergyD5(const std::vector<double> &lambda_pows, const std::vector<double> &dist_pows);
template double XCEnergyD6(const std::vector<double> &lambda_pows, const std::vector<double> &dist_pows);
template double XCEnergyD7(const std::vector<double> &lambda_pows, const std::vector<double> &dist_pows);
template double XCEnergyD8(const std::vector<double> &lambda_pows, const std::vector<double> &dist_pows);
template double XCEnergyD9(const std::vector<double> &lambda_pows, const std::vector<double> &dist_pows);
template double XCEnergyD10(const std::vector<double> &lambda_pows, const std::vector<double> &dist_pows);

template autodiff::dual XCEnergy1(const std::vector<double> &lambda_pows, const std::vector<autodiff::dual> &dist_pows);
template autodiff::dual XCEnergy2(const std::vector<double> &lambda_pows, const std::vector<autodiff::dual> &dist_pows);
template autodiff::dual XCEnergy3(const std::vector<double> &lambda_pows, const std::vector<autodiff::dual> &dist_pows);
template autodiff::dual XCEnergy4(const std::vector<double> &lambda_pows, const std::vector<autodiff::dual> &dist_pows);
template autodiff::dual XCEnergy5(const std::vector<double> &lambda_pows, const std::vector<autodiff::dual> &dist_pows);
template autodiff::dual XCEnergy6(const std::vector<double> &lambda_pows, const std::vector<autodiff::dual> &dist_pows);
template autodiff::dual XCEnergy7(const std::vector<double> &lambda_pows, const std::vector<autodiff::dual> &dist_pows);
template autodiff::dual XCEnergy8(const std::vector<double> &lambda_pows, const std::vector<autodiff::dual> &dist_pows);
template autodiff::dual XCEnergy9(const std::vector<double> &lambda_pows, const std::vector<autodiff::dual> &dist_pows);
template autodiff::dual XCEnergy10(const std::vector<double> &lambda_pows, const std::vector<autodiff::dual> &dist_pows);

template autodiff::dual XCEnergyD1(const std::vector<double> &lambda_pows, const std::vector<autodiff::dual> &dist_pows);
template autodiff::dual XCEnergyD2(const std::vector<double> &lambda_pows, const std::vector<autodiff::dual> &dist_pows);
template autodiff::dual XCEnergyD3(const std::vector<double> &lambda_pows, const std::vector<autodiff::dual> &dist_pows);
template autodiff::dual XCEnergyD4(const std::vector<double> &lambda_pows, const std::vector<autodiff::dual> &dist_pows);
template autodiff::dual XCEnergyD5(const std::vector<double> &lambda_pows, const std::vector<autodiff::dual> &dist_pows);
template autodiff::dual XCEnergyD6(const std::vector<double> &lambda_pows, const std::vector<autodiff::dual> &dist_pows);
template autodiff::dual XCEnergyD7(const std::vector<double> &lambda_pows, const std::vector<autodiff::dual> &dist_pows);
template autodiff::dual XCEnergyD8(const std::vector<double> &lambda_pows, const std::vector<autodiff::dual> &dist_pows);
template autodiff::dual XCEnergyD9(const std::vector<double> &lambda_pows, const std::vector<autodiff::dual> &dist_pows);
template autodiff::dual XCEnergyD10(const std::vector<double> &lambda_pows, const std::vector<autodiff::dual> &dist_pows);

// Explicit instantiation of all the types for the naive model functional.
template double MolFFSim::NaiveModelE(const double lambda, const double &dist);
template autodiff::dual MolFFSim::NaiveModelE(const double lambda,const autodiff::dual &dist);

template double MolFFSim::XCSpheSymm(const double lambda, const double &dist);
template autodiff::dual MolFFSim::XCSpheSymm(const double lambda,const autodiff::dual &dist);

template double MolFFSim::XCCylinSymm(const double lambda, const double &dist);
template autodiff::dual MolFFSim::XCCylinSymm(const double lambda,const autodiff::dual &dist);
