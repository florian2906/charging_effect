#include <iostream>
#include <filesystem> 
#include <viennaps.hpp>
#include <geometries/psMakeTrench.hpp>
#include <random>

#include "Ion.hpp"
#include "GeometryCollider.hpp"
#include "PoissonSolver.hpp"


int main(int argc, char** argv) {
    std::cout << "\n==================================================" << std::endl;
    std::cout << "       VIENNAPS CHARGING EFFECT SIMULATION        " << std::endl;
    std::cout << "==================================================\n" << std::endl;
    



    std::cout << "[1/4] INITIALIZATION & CONFIGURATION" << std::endl;
    std::cout << "--------------------------------------------------" << std::endl;
    constexpr int DIM = 2;

    // Create the output directory if it does not already exist
    if (!std::filesystem::exists("output")) {
        std::cout << "  [*] Creating output directory 'output/'..." << std::endl;
        std::filesystem::create_directory("output");
    }

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



    
    std::cout << "[2/4] GEOMETRY GENERATION" << std::endl;
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

    // Initialize grid memory and neighborhood immediately
    std::cout << "  [*] Building cell neighborhood and allocating fields..." << std::endl;
    auto cellSet = domain->getCellSet();
    cellSet->buildNeighborhood();
    cellSet->addScalarData("ChargeDensity", 0.0);
    cellSet->addScalarData("ElectricPotential", 0.0);
    cellSet->addScalarData("E_x", 0.0);
    cellSet->addScalarData("E_y", 0.0);

    // Export Surface Mesh for visualization
    std::cout << "  [*] Extracting and exporting surface mesh..." << std::endl;
    auto surfaceMesh = domain->getDiskMesh(); 
    viennals::VTKWriter<double>(surfaceMesh, "01_trench_surface.vtp").apply();
    std::cout << "      -> Saved: '01_trench_surface.vtp'\n" << std::endl;




    std::cout << "[3/4] MONTE CARLO CHARGE DEPOSITION" << std::endl;
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
    std::vector<Ion> all_ions;

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
        all_ions.push_back(ion);
    }

    std::cout << "      -> " << hit_count << " ions successfully hit the surface." << std::endl;

    // Export the charge distribution
    std::cout << "  [*] Exporting charge distribution to ParaView..." << std::endl;
    surfaceMesh->getCellData().insertNextScalarData(surfaceCharge, "SurfaceCharge");
    viennals::VTKWriter<double>(surfaceMesh, "02_trench_with_charge.vtp").apply();
    std::cout << "      -> Saved: '02_trench_with_charge.vtp'\n" << std::endl;

    std::cout << "  [*] Exporting sample ion trajectories..." << std::endl;
    auto pathMesh = viennals::Mesh<double>::New();
    
    // Export the first 50 ions to visualize their paths
    int paths_to_export = std::min(num_ions, 50); 
    
    for (int i = 0; i < paths_to_export; ++i) {
        const auto& ion_path = all_ions[i].trajectory; 
        
        unsigned int current_nodes = pathMesh->nodes.size();
        
        for (int j = 0; j < ion_path.size(); ++j) {
            pathMesh->nodes.push_back(ion_path[j]);
            if (j > 0) {
                pathMesh->insertNextLine({current_nodes + j - 1, current_nodes + j});
            }
        }
    }
    viennals::VTKWriter<double>(pathMesh, "03_ion_trajectories.vtp").apply();
    std::cout << "      -> Saved: '03_ion_trajectories.vtp'\n" << std::endl;




    std::cout << "\n[4/4] POISSON SOLVER & FIELD MAPPING" << std::endl;
    std::cout << "--------------------------------------------------" << std::endl;

    std::cout << "  [*] Mapping surface charge to volumetric CellSet..." << std::endl;

    // Get the total number of cells
    int num_cells = cellSet->getElements().size();

    // Just GET the memory pointer, it is already perfectly sized!
    auto chargeDensityPtr = cellSet->getScalarData("ChargeDensity");
    
    // Extract the surface normals to push the charge into the gas
    auto normals = surfaceMesh->getCellData().getVectorData("Normals");

    double grid_delta = params.get("gridDelta");
    double cell_volume = grid_delta * grid_delta; 

    // Loop over all surface nodes to map their charge into the 2D grid
    for (size_t i = 0; i < surfaceCharge.size(); ++i) {
        if (surfaceCharge[i] > 0.0) {
            
            std::array<double, 3> node_coord = surfaceMesh->nodes[i];
            std::array<double, 3> normal = normals->at(i);

            // Push the coordinate half a grid-cell outward into the gas
            node_coord[0] += normal[0] * (grid_delta * 0.5);
            node_coord[1] += normal[1] * (grid_delta * 0.5);

            int cell_idx = cellSet->getIndex(node_coord);
            
            // If the cell exists, deposit the charge density directly into ViennaPS's memory
            if (cell_idx >= 0) {
                // Use (*chargeDensityPtr) to access the underlying vector
                (*chargeDensityPtr)[cell_idx] += surfaceCharge[i] / cell_volume;
            }
        }
    }

    // Call solver to compute potential and fields
    std::cout << "  [*] Solving Poisson Equation for Electric Potential..." << std::endl;
    PoissonSolver::solvePoisson(domain, grid_delta);

    std::cout << "  [*] Exporting volumetric grid with mapped fields..." << std::endl;
    domain->getCellSet()->writeVTU("04_volume_grid_for_poisson.vtu");
    std::cout << "      -> Saved: '04_volume_grid_for_poisson.vtu'\n" << std::endl;


    std::cout << "==================================================" << std::endl;
    std::cout << "               SIMULATION COMPLETE                " << std::endl;
    std::cout << "==================================================\n" << std::endl;
    return 0;
}