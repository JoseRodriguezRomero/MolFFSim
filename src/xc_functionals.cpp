#define M_PI_SQRT                   1.7724538509055159
#define MIN_ATOM_POTENTIAL_DIST     1.0E-6

#define POW_TABLE_SIZE              18

// Solve for x "log(1-erf(sqrt(x)))/log(10) = -12"
#define APPROX_RADIUS               25.4220809376487

#include "xc_functionals.hpp"

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
        // If this is true then the naive potential part can be approximated
        // as 1/r. Therefore, its Laplacian must be approximately zero.
        return T(0.0);
    }
    
    std::vector<T> dist_pows;
    std::vector<double> lambda_pows;
    
    dist_pows.reserve(POW_TABLE_SIZE);
    lambda_pows.reserve(POW_TABLE_SIZE);
    
    dist_pows[0] = dist;
    lambda_pows[0] = lambda;
    
    if (dist == 0) {
        for (unsigned i = 1; i < POW_TABLE_SIZE; i++) {
            dist_pows[i] = 0;
            lambda_pows[i] = lambda_pows[i-1] * lambda;
        }
    }
    else {
        for (unsigned i = 1; i < POW_TABLE_SIZE; i++) {
            dist_pows[i] = dist_pows[i-1] * dist;
            lambda_pows[i] = lambda_pows[i-1] * lambda;
        }
    }
    
    // THIS IS PROBLEMATIC WITH AUTODIFF!!!
    // T aux_c = sqrt(lambda) /  (exp(lambda*dist_pows[2]) * M_PI_SQRT);
    T aux_c = sqrt(lambda) * exp(-lambda*dist_pows[1]) / M_PI_SQRT;
    
    // This expression was fully factored using Maxima, and obtained from an
    // order ten XC expansion (in principle this preferable performance-wise).
    T xc_energy = -((lambda_pows[0]*(512*dist_pows[17]*lambda_pows[17]-
        43776*dist_pows[15]*lambda_pows[16]-48640*dist_pows[15]*lambda_pows[15]+
        1488384*dist_pows[13]*lambda_pows[15]+3307520*dist_pows[13]*
        lambda_pows[14]-26046720*dist_pows[11]*lambda_pows[14]+3720960*
        dist_pows[13]*lambda_pows[13]-86822400*dist_pows[11]*lambda_pows[13]+
        253955520*dist_pows[9]*lambda_pows[13]-195350400*dist_pows[11]*
        lambda_pows[12]+1128691200*dist_pows[9]*lambda_pows[12]-1396755360*
        dist_pows[7]*lambda_pows[12]-223257600*dist_pows[11]*lambda_pows[11]+
        3809332800*dist_pows[9]*lambda_pows[11]-7759752000*dist_pows[7]*
        lambda_pows[11]+4190266080*dist_pows[5]*lambda_pows[11]+8707046400*
        dist_pows[9]*lambda_pows[10]-34918884000*dist_pows[7]*lambda_pows[10]+
        27935107200*dist_pows[5]*lambda_pows[10]-6285399120*dist_pows[3]*
        lambda_pows[10]+10158220800*dist_pows[9]*lambda_pows[9]-119721888000*
        dist_pows[7]*lambda_pows[9]+157134978000*dist_pows[5]*lambda_pows[9]-
        48886437600*dist_pows[3]*lambda_pows[9]+3928374450*dist_pows[1]*
        lambda_pows[9]-279351072000*dist_pows[7]*lambda_pows[8]+718331328000*
        dist_pows[5]*lambda_pows[8]-329983453800*dist_pows[3]*lambda_pows[8]+
        34918884000*dist_pows[1]*lambda_pows[8]-654729075*lambda_pows[8]-
        335221286400*dist_pows[7]*lambda_pows[7]+2514159648000*dist_pows[5]*
        lambda_pows[7]-1885619736000*dist_pows[3]*lambda_pows[7]+274986211500*
        dist_pows[1]*lambda_pows[7]-6547290750*lambda_pows[7]+6033983155200*
        dist_pows[5]*lambda_pows[6]-8799558768000*dist_pows[3]*lambda_pows[6]+
        1885619736000*dist_pows[1]*lambda_pows[6]-58925616750*lambda_pows[6]+
        7542478944000*dist_pows[5]*lambda_pows[5]-31678411564800*dist_pows[3]*
        lambda_pows[5]+10999448460000*dist_pows[1]*lambda_pows[5]-471404934000*
        lambda_pows[5]-79196028912000*dist_pows[3]*lambda_pows[4]+
        52797352608000*dist_pows[1]*lambda_pows[4]-3299834538000*lambda_pows[4]-
        105594705216000*dist_pows[3]*lambda_pows[3]+197990072280000*
        dist_pows[1]*lambda_pows[3]-19799007228000*lambda_pows[3]+
        527973526080000*dist_pows[1]*lambda_pows[2]-98995036140000*
        lambda_pows[2]+791960289120000*dist_pows[1]*lambda_pows[1]-
        395980144560000*lambda_pows[1]-1187940433680000*lambda_pows[0]-
        2375880867360000))/1187940433680000);
        
    return aux_c * xc_energy;
}

template<typename T>
T MolFFSim::XCCylinSymm(const double lambda, const T &dist) {
    if ((lambda*dist*dist) >= APPROX_RADIUS) {
        // If this is true then the naive potential part can be approximated
        // as 1/r. Therefore, its Laplacian must be approximately zero.
        return T(0.0);
    }
    
    std::vector<T> dist_pows;
    std::vector<double> lambda_pows;
    
    dist_pows.reserve(POW_TABLE_SIZE);
    lambda_pows.reserve(POW_TABLE_SIZE);
    
    dist_pows[0] = dist;
    lambda_pows[0] = lambda;
    
    if (dist == 0) {
        for (unsigned i = 1; i < POW_TABLE_SIZE; i++) {
            dist_pows[i] = 0;
            lambda_pows[i] = lambda_pows[i-1] * lambda;
        }
    }
    else {
        for (unsigned i = 1; i < POW_TABLE_SIZE; i++) {
            dist_pows[i] = dist_pows[i-1] * dist;
            lambda_pows[i] = lambda_pows[i-1] * lambda;
        }
    }
    
    // THIS IS PROBLEMATIC WITH AUTODIFF!!!
    // T aux_c = sqrt(lambda) /  (exp(lambda*dist_pows[2]) * M_PI_SQRT);
    T aux_c = sqrt(lambda) * exp(-lambda*dist_pows[1]) / M_PI_SQRT;
    
    // This expression was fully factored using Maxima, and obtained from an
    // order ten XC expansion (in principle this preferable performance-wise).
    T xc_energy = -((dist_pows[0]*lambda_pows[1]*(512*dist_pows[17]*
        lambda_pows[17]-48384*dist_pows[15]*lambda_pows[16]-53760*dist_pows[15]*
        lambda_pows[15]+1838592*dist_pows[13]*lambda_pows[15]+4085760*
        dist_pows[13]*lambda_pows[14]-36465408*dist_pows[11]*lambda_pows[14]+
        4596480*dist_pows[13]*lambda_pows[13]-121551360*dist_pows[11]*
        lambda_pows[13]+410235840*dist_pows[9]*lambda_pows[13]-273490560*
        dist_pows[11]*lambda_pows[12]+1823270400*dist_pows[9]*lambda_pows[12]-
        2666532960*dist_pows[7]*lambda_pows[12]-312560640*dist_pows[11]*
        lambda_pows[11]+6153537600*dist_pows[9]*lambda_pows[11]-14814072000*
        dist_pows[7]*lambda_pows[11]+9777287520*dist_pows[5]*lambda_pows[11]+
        14065228800*dist_pows[9]*lambda_pows[10]-66663324000*dist_pows[7]*
        lambda_pows[10]+65181916800*dist_pows[5]*lambda_pows[10]-18856197360*
        dist_pows[3]*lambda_pows[10]+16409433600*dist_pows[9]*lambda_pows[9]-
        228559968000*dist_pows[7]*lambda_pows[9]+366648282000*dist_pows[5]*
        lambda_pows[9]-146659312800*dist_pows[3]*lambda_pows[9]+16499172690*
        dist_pows[1]*lambda_pows[9]-533306592000*dist_pows[7]*lambda_pows[8]+
        1676106432000*dist_pows[5]*lambda_pows[8]-989950361400*dist_pows[3]*
        lambda_pows[8]+146659312800*dist_pows[1]*lambda_pows[8]-4583103525*
        lambda_pows[8]-639967910400*dist_pows[7]*lambda_pows[7]+5866372512000*
        dist_pows[5]*lambda_pows[7]-5656859208000*dist_pows[3]*lambda_pows[7]+
        1154942088300*dist_pows[1]*lambda_pows[7]-45831035250*lambda_pows[7]+
        14079294028800*dist_pows[5]*lambda_pows[6]-26398676304000*dist_pows[3]*
        lambda_pows[6]+7919602891200*dist_pows[1]*lambda_pows[6]-412479317250*
        lambda_pows[6]+17599117536000*dist_pows[5]*lambda_pows[5]-
        95035234694400*dist_pows[3]*lambda_pows[5]+46197683532000*dist_pows[1]*
        lambda_pows[5]-3299834538000*lambda_pows[5]-237588086736000*
        dist_pows[3]*lambda_pows[4]+221748880953600*dist_pows[1]*lambda_pows[4]-
        23098841766000*lambda_pows[4]-316784115648000*dist_pows[3]*
        lambda_pows[3]+831558303576000*dist_pows[1]*lambda_pows[3]-
        138593050596000*lambda_pows[3]+2217488809536000*dist_pows[1]*
        lambda_pows[2]-692965252980000*lambda_pows[2]+3326233214304000*
        dist_pows[1]*lambda_pows[1]-2771861011920000*lambda_pows[1]-
        8315583035760000*lambda_pows[0]-16631166071520000))/12473374553640000);
    
    return aux_c * xc_energy;
}

// Explicit instantiation of all the types for the naive model functional.
template double MolFFSim::NaiveModelE(const double lambda, const double &dist);
template autodiff::dual MolFFSim::NaiveModelE(const double lambda,const autodiff::dual &dist);

template double MolFFSim::XCSpheSymm(const double lambda, const double &dist);
template autodiff::dual MolFFSim::XCSpheSymm(const double lambda,const autodiff::dual &dist);

template double MolFFSim::XCCylinSymm(const double lambda, const double &dist);
template autodiff::dual MolFFSim::XCCylinSymm(const double lambda,const autodiff::dual &dist);
