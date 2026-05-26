#include <iostream>
#include <filesystem> // Required for changing directories
#include <viennaps.hpp>
#include <geometries/psMakeTrench.hpp>
#include <random>

#include "Ion.hpp"
#include "GeometryCollider.hpp"


int main(int argc, char** argv) {
    std::cout << "\n==================================================" << std::endl;
    std::cout << "       VIENNAPS CHARGING EFFECT SIMULATION        " << std::endl;
    std::cout << "==================================================\n" << std::endl;
    
    std::cout << "[1/3] INITIALIZATION & CONFIGURATION" << std::endl;
    std::cout << "--------------------------------------------------" << std::endl;

    // Create the output directory if it does not already exist
    if (!std::filesystem::exists("output")) {
        std::cout << "  [*] Creating output directory 'output/'..." << std::endl;
        std::filesystem::create_directory("output");
    }

    constexpr int DIM = 2;

    // Load Parameters from configuration file
    viennaps::util::Parameters params;
    try {
        params.readConfigFile("Parameters.txt");
        std::cout << "  [*] Simulation parameters loaded:" << std::endl;
        for (const auto& [key, value] : params.m) {
            std::cout << "      - " << key << ": " << value << std::endl;}
    }
    catch (const std::exception& e) {
        std::cerr << "  [!] Error loading parameters: " << e.what() << std::endl;
        return 1;
    }

    // Change working directory to output for all following file operations
    std::filesystem::current_path("output");
    std::cout << "  [+] Working directory changed to ./output/\n" << std::endl;

    
    std::cout << "[2/3] GEOMETRY GENERATION" << std::endl;
    std::cout << "--------------------------------------------------" << std::endl;

    // Create Domain and Trench
    std::cout << "  [*] Setting up the simulation domain..." << std::endl;
    auto domain = viennaps::Domain<double, DIM>::New(
        params.get("gridDelta"),    
        params.get("xExtent"),      
        viennaps::BoundaryType::REFLECTIVE_BOUNDARY
    );

    std::cout << "  [*] Applying trench geometry..." << std::endl;
    viennaps::MakeTrench<double, DIM>(
        domain, 
        params.get("trenchWidth"),                
        params.get("trenchHeight"), 
        params.get("taperAngle"),
        params.get("maskHeight")
    ).apply();

    // Generate and Export Cell Set
    std::cout << "  [*] Generating cell set..." << std::endl;
    domain->generateCellSet(
        params.get("maskHeight") + 10 * params.get("gridDelta"), 
        viennaps::Material::GAS, 
        true    
    );

    // Export Surface Mesh for visualization
    std::cout << "  [*] Extracting and exporting surface mesh..." << std::endl;
    auto surfaceMesh = domain->getDiskMesh(); 
    viennals::VTKWriter<double>(surfaceMesh, "01_trench_surface.vtp").apply();
    std::cout << "      -> Saved: '01_trench_surface.vtp'" << std::endl;

    std::cout << "  [*] Exporting volumetric grid for Poisson solver..." << std::endl;
    domain->getCellSet()->writeVTU("02_volume_grid_for_poisson.vtu");
    std::cout << "      -> Saved: '02_volume_grid_for_poisson.vtu'\n" << std::endl;


    std::cout << "[3/3] MONTE CARLO CHARGE DEPOSITION" << std::endl;
    std::cout << "--------------------------------------------------" << std::endl;

    // Set up geometry collider
    std::cout << "  [*] Building Embree collision geometry..." << std::endl;
    GeometryCollider<DIM> collider;
    collider.buildSurface(surfaceMesh, params.get("gridDelta"));

    std::vector<double> surfaceCharge(surfaceMesh->nodes.size(), 0.0);
    int num_ions = params.get<int>("numRaysPerPoint") * 10;

    double dt = 0.01;
    double ion_mass = params.get("mass");
    double ion_maxVel = params.get("maxVel");
    double x_extent = params.get("xExtent");
    double start_y = params.get("maskHeight") + 10.0;
    double y_kill_bound = -(params.get("trenchHeight")) - 5.0;
    double x_kill_bound = x_extent / 2.0;

    std::cout << "  [*] Launching a batch of " << num_ions << " ions..." << std::endl;
    int hit_count = 0;

    std::mt19937 gen(std::random_device{}());
    std::uniform_real_distribution<double> dist_x(-(x_extent / 2.0), (x_extent / 2.0));

    for (int i = 0; i < num_ions; ++i) {
        // Initialize an Ion at the top of the domain at a random X position, pointing down
        double start_x = dist_x(gen);
        Ion ion(ion_mass, 1.0, {start_x, start_y, 0.0}, {0.0, -ion_maxVel, 0.0});

        bool active = true;

        // Ray tracing loop
        while (active) {
            std::array<double, 3> old_pos = ion.position;

            double E_x = 0.0;
            double E_y = 0.0;

            ion.advanceTimeStep(dt, E_x, E_y);

            double dx = ion.position[0] - old_pos[0];
            double dy = ion.position[1] - old_pos[1];
            double step_distance = sqrt(dx*dx + dy*dy);

            std::array<double, 3> direction = {dx/step_distance, dy/step_distance, 0.0};
        
            int hit_id = collider.checkCollision(old_pos, direction, step_distance);

            // Handle Collision
            if (hit_id >= 0) {
                surfaceCharge[hit_id] += ion.charge; // Deposit Charge
                hit_count++;
                active = false; // Kill the particle
            }
            
            // Bounds checking (Kill particle if it flies out of the domain)
            if (ion.position[1] < y_kill_bound || std::abs(ion.position[0]) > x_kill_bound) {
                active = false; 
            }
        }
    }

    std::cout << "      -> " << hit_count << " ions successfully hit the surface." << std::endl;

    // Export the charge distribution
    std::cout << "  [*] Exporting charge distribution to ParaView..." << std::endl;
    surfaceMesh->getCellData().insertNextScalarData(surfaceCharge, "SurfaceCharge");
    viennals::VTKWriter<double>(surfaceMesh, "03_trench_with_charge.vtp").apply();
    std::cout << "      -> Saved: '03_trench_with_charge.vtp'\n" << std::endl;


    std::cout << "==================================================" << std::endl;
    std::cout << "               SIMULATION COMPLETE                " << std::endl;
    std::cout << "==================================================\n" << std::endl;

    return 0;
}