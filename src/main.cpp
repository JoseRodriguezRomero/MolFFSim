#define MAX_PRINT_BUFFER_SIZE       256

#include <fstream>
#include <iostream>

#include <autodiff/reverse/var.hpp>
#include <autodiff/forward/dual.hpp>

#include "system.hpp"

int main(int argc, char** argv) {
    // True if forward-mode automatic differentiation is to be used.
    bool diff_is_forward = false;
    
    // True if a single-point energy calculation is to be made.
    bool point_energy = true;
    
    // True if the calculation consists of displacing and rotating all the
    // molecules in the system until an energy minimum is reached.
    bool system_geom_optim = false;
    
    // True if the calculation consists of adjusting the internal coordinates
    // of all the atoms of each molecule until their monomer energy reach a
    // minimum of energy. The data provided in the MOLECULE_LIST section of
    // any input file is unused in the mode of calculation.
    bool molec_geom_optim = false;
    
    if (argc < 2) {
        std::cout << "Please specify an input file." << std::endl;
        return 0;
    }
    
    std::ifstream input_file;
    input_file.open(argv[1], std::fstream::in);
    
    for (unsigned i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "point_energy")) {
            point_energy = true;
            system_geom_optim = false;
            molec_geom_optim = false;
        }
        else if (!strcmp(argv[i], "system_geom_optim")) {
            point_energy = false;
            system_geom_optim = true;
            molec_geom_optim = false;
        }
        else if (!strcmp(argv[i], "molec_geom_optim")) {
            point_energy = false;
            system_geom_optim = false;
            molec_geom_optim = true;
        }
        else if (!strcmp(argv[i], "forward_mode")) {
            diff_is_forward = true;
        }
        else if (!strcmp(argv[i], "reverse_mode")) {
            diff_is_forward = false;
        }
    }
    
    char buffer[MAX_PRINT_BUFFER_SIZE];
    
    if (point_energy) {
        MolFFSim::System<double> system;
        system.ReadInputFile(input_file);
        input_file.close();
        
        std::cout << system << std::endl;;
        
        snprintf(buffer,MAX_PRINT_BUFFER_SIZE, "%-20s %-15.5E [kJ/mol]",
                 "Total Energy:",0.0);
        std::cout << buffer << std::endl;
        
        snprintf(buffer,MAX_PRINT_BUFFER_SIZE, "%-20s %-15.5E [kJ/mol]",
                 "Interaction Energy:",0.0);
        std::cout << buffer << std::endl << std::endl;
    }
    else if (diff_is_forward) {
        MolFFSim::System<autodiff::dual> system;
        system.ReadInputFile(input_file);
        input_file.close();
        
        if (system_geom_optim) {
            std::cout << system << std::endl;;
            
            snprintf(buffer,MAX_PRINT_BUFFER_SIZE, "%-20s %-15.5E [kJ/mol]",
                     "Total Energy:",0.0);
            std::cout << buffer << std::endl;
            
            snprintf(buffer,MAX_PRINT_BUFFER_SIZE, "%-20s %-15.5E [kJ/mol]",
                     "Interaction Energy:",0.0);
            std::cout << buffer << std::endl << std::endl;
        }
        else if (molec_geom_optim) {
        }
    }
    else {
        MolFFSim::System<autodiff::var> system;
        system.ReadInputFile(input_file);
        input_file.close();
        
        if (system_geom_optim) {
            std::cout << system << std::endl;;
            
            snprintf(buffer,MAX_PRINT_BUFFER_SIZE, "%-20s %-15.5E [kJ/mol]",
                     "Total Energy:",double(system.SystemEnergy()));
            std::cout << buffer << std::endl;
            
            snprintf(buffer,MAX_PRINT_BUFFER_SIZE, "%-20s %-15.5E [kJ/mol]",
                     "Interaction Energy:",
                     double(system.SystemInteractionEnergy()));
            std::cout << buffer << std::endl << std::endl;
        }
        else if (molec_geom_optim) {
        }
    }
    
    return 0;
}
