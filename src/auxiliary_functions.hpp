#ifndef AUX_FUNCTIONS
#define AUX_FUNCTIONS

#include <string>

namespace MolFFSim {

std::string LabelFromAtomicNumber(const unsigned atomic_number);
unsigned AtomicNumberFromLabel(const std::string &element);
unsigned ECPEffectiveAtomicNumber(const unsigned atomic_number);

}

#endif
