#include <iostream>
#include <filesystem> // Required for changing directories
#include <viennaps.hpp>
#include <geometries/psMakeTrench.hpp>


int main(int argc, char** argv) {
    std::cout << "--- Initializing Charging Simulation ---" << " \n " << std::endl;
    
    // Create the output directory
    std::filesystem::create_directory("output");

    constexpr int DIM = 2;

    // Load Parameters from configuration file
    viennaps::util::Parameters params;
    try {
        params.readConfigFile("Parameters.txt");
        std::cout << "Simulation parameters loaded:" << std::endl;
        for (const auto& [key, value] : params.m) {
            std::cout << "  " << key << ": " << value << std::endl;}
    }
    catch (const std::exception& e) {
        std::cerr << "Error loading parameters: " << e.what() << std::endl;
        return 1;
    }
    std::cout << "-----------------------------------------" << " \n " << std::endl;

    // Change working directory to output for all following file operations
    std::filesystem::current_path("output");
    std::cout << "Working directory changed to ./output/\n" << std::endl;

    // Create Domain and Trench
    std::cout << "Setting up the simulation domain..." << std::endl;
    auto domain = viennaps::Domain<double, DIM>::New(
        params.get("gridDelta"),    
        params.get("xExtent"),      
        viennaps::BoundaryType::REFLECTIVE_BOUNDARY
    );

    std::cout << "Applying trench geometry..." << std::endl;
    viennaps::MakeTrench<double, DIM>(
        domain, 
        params.get("trenchWidth"),                
        params.get("trenchHeight"), 
        params.get("taperAngle"),
        params.get("maskHeight")
    ).apply();

    // Export Surface Mesh for visualization
    std::cout << "Extracting and exporting surface mesh..." << std::endl;
    auto surfaceMesh = domain->getDiskMesh(); 
    viennals::VTKWriter<double>(surfaceMesh, "01_trench_surface.vtp").apply();

    // Generate and Export Cell Set
    std::cout << "Generating cell set..." << std::endl;
    domain->generateCellSet(
        params.get("maskHeight") + 10 * params.get("gridDelta"), 
        viennaps::Material::GAS, 
        true    
    );
    std::cout << "Exporting volumetric grid for Poisson solver..." << std::endl;
    domain->getCellSet()->writeVTU("02_volume_grid_for_poisson.vtu");
    
    std::cout << "-----------------------------------------" << " \n " << std::endl;

    return 0;
}