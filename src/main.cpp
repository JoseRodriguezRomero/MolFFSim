#define MAX_PRINT_BUFFER_SIZE       256
#define HARTREE_TO_KJ_MOL           2625.5002

#include <fstream>
#include <iostream>
#include <filesystem>

#include "omp.h"

#include <Eigen/Core>
#include <autodiff/forward/dual.hpp>

#include "system.hpp"

int main(int argc, char** argv) {
    int n = std::thread::hardware_concurrency();
    omp_set_num_threads(n);
    Eigen::setNbThreads(n);
    
    Eigen::initParallel();
    
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
    
    // True if the forces of each atom are to be printed.
    bool point_energy_forces = false;
    
    if (argc < 2) {
        std::cout << "Please specify an input file." << std::endl;
        return 0;
    }
    
    std::ifstream input_file;
    input_file.open(argv[1], std::fstream::in);
    
    std::string backup_filename1 = "backup1_";
    backup_filename1 += std::filesystem::path(argv[1]).filename();
    backup_filename1 =
        std::filesystem::path(argv[1]).replace_filename(backup_filename1);
    
    std::string backup_filename2 = "backup2_";
    backup_filename2 += std::filesystem::path(argv[1]).filename();
    backup_filename2 =
        std::filesystem::path(argv[1]).replace_filename(backup_filename2);
    
    for (unsigned i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "point_energy")) {
            point_energy = true;
            point_energy_forces = false;
            system_geom_optim = false;
            molec_geom_optim = false;
        }
        if (!strcmp(argv[i], "point_energy_forces")) {
            point_energy = false;
            point_energy_forces = true;
            system_geom_optim = false;
            molec_geom_optim = false;
        }
        else if (!strcmp(argv[i], "system_geom_optim")) {
            point_energy = false;
            point_energy_forces = false;
            system_geom_optim = true;
            molec_geom_optim = false;
        }
        else if (!strcmp(argv[i], "molec_geom_optim")) {
            point_energy = false;
            point_energy_forces = false;
            system_geom_optim = false;
            molec_geom_optim = true;
        }
    }
    
    char buffer[MAX_PRINT_BUFFER_SIZE];
    
    if (point_energy) {
        MolFFSim::System<double> system;
        system.ReadInputFile(input_file);
        input_file.close();
        
        system.PolarizeMolecules();
        std::cout << system << std::endl;;
                
        snprintf(buffer,MAX_PRINT_BUFFER_SIZE, "%-20s %15.5E [kJ/mol]",
                 "Total Energy:", HARTREE_TO_KJ_MOL * system.SystemEnergy());
        std::cout << buffer << std::endl;
        
        snprintf(buffer,MAX_PRINT_BUFFER_SIZE, "%-20s %15.5E [kJ/mol]",
                 "Interaction Energy:", 
                 HARTREE_TO_KJ_MOL * system.SystemInteractionEnergy());
        std::cout << buffer << std::endl << std::endl;
    }
    else if (point_energy_forces) {
        MolFFSim::System<autodiff::dual> system;
        system.ReadInputFile(input_file);
        input_file.close();
        
        system.PolarizeMolecules();
        std::cout << system << std::endl;;
                
        snprintf(buffer,MAX_PRINT_BUFFER_SIZE, "%-20s %15.5E [kJ/mol]",
                 "Total Energy:", HARTREE_TO_KJ_MOL *
                 double(system.SystemEnergy()));
        std::cout << buffer << std::endl;
        
        snprintf(buffer,MAX_PRINT_BUFFER_SIZE, "%-20s %15.5E [kJ/mol]",
                 "Interaction Energy:", HARTREE_TO_KJ_MOL *
                 double(system.SystemInteractionEnergy()));
        std::cout << buffer << std::endl << std::endl;
    }
    else {
        MolFFSim::System<autodiff::dual> system;
        system.ReadInputFile(input_file);
        input_file.close();
        
        system.setBackupFile1(backup_filename1);
        system.setBackupFile2(backup_filename2);
        
        if (system_geom_optim) {
            system.PolarizeMolecules();
            system.OptimizeGeometry(std::cout);
            
            std::cout << system << std::endl;;
            
            snprintf(buffer,MAX_PRINT_BUFFER_SIZE, "%-20s %15.5E [kJ/mol]",
                     "Total Energy:",
                     HARTREE_TO_KJ_MOL * double(system.SystemEnergy()));
            std::cout << buffer << std::endl;
            
            snprintf(buffer,MAX_PRINT_BUFFER_SIZE, "%-20s %15.5E [kJ/mol]",
                     "Interaction Energy:", HARTREE_TO_KJ_MOL *
                     double(system.SystemInteractionEnergy()));
            std::cout << buffer << std::endl << std::endl;
        }
        else if (molec_geom_optim) {
            system.OptimizeMoleculeGeometries();
        }
    }
    
    return 0;
}
