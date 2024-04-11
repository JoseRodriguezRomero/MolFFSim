#include "auxiliary_functions.hpp"

std::string MolFFSim::LabelFromAtomicNumber(const unsigned atomic_number) {
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

unsigned MolFFSim::AtomicNumberFromLabel(const std::string &element) {
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
    
    // Unsupported elements have atomic number equal to zero.
    return 0;
}

unsigned MolFFSim::ECPEffectiveAtomicNumber(const unsigned atomic_number) {
    switch (atomic_number) {
        case 1:
            return 1;
            break;
        case 2:
            return 2;
            break;
        case 3:
            return 1;
            break;
        case 4:
            return 2;
            break;
        case 5:
            return 3;
            break;
        case 6:
            return 4;
            break;
        case 7:
            return 5;
            break;
        case 8:
            return 6;
            break;
        case 9:
            return 7;
            break;
        case 10:
            return 8;
            break;
        case 11:
            return 1;
            break;
        case 12:
            return 2;
            break;
        case 13:
            return 3;
            break;
        case 14:
            return 4;
            break;
        case 15:
            return 5;
            break;
        case 16:
            return 6;
            break;
        case 17:
            return 7;
            break;
        case 18:
            return 8;
            break;
        default:
            break;
    }
    
    // Unsupported elements have atomic number equal to zero.
    return 0;
}
