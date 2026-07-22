#include <iostream>
#include <filesystem> 
#include <viennaps.hpp>
#include <geometries/psMakeTrench.hpp>
#include <random>

#include "Particle.hpp"
#include "GeometryCollider.hpp"
#include "PoissonSolver.hpp"
#include "Helpers.hpp"
#include "EtchVelocityField.hpp"


int main(int argc, char** argv) {
    std::cout << "\n========================================================================" << std::endl;
    std::cout << "                  VIENNAPS CHARGING EFFECT SIMULATION                   " << std::endl;
    std::cout << "========================================================================\n" << std::endl;
    


    // ========================================================================
    // [PHASE 1] INITIALIZATION & CONFIGURATION
    // ========================================================================
    std::cout << "[PHASE 1] INITIALIZATION & CONFIGURATION" << std::endl;
    std::cout << "------------------------------------------------------------------------" << std::endl;
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
        std::cout << "  [*] Simulation parameters loaded successfully:" << std::endl;
        for (const auto& [key, value] : params.m) {
            std::cout << "      - " << key << ": " << value << std::endl;
        }
    }
    catch (const std::exception& e) {
        std::cerr << "  [!] Error loading parameters: " << e.what() << std::endl;
        return 1;
    }

    // Change working directory to output for all following file operations
    std::filesystem::current_path("output");
    std::cout << "  [+] Working directory changed to ./output/\n" << std::endl;



    // ========================================================================
    // [PHASE 2] INITIAL GEOMETRY GENERATION
    // ========================================================================
    std::cout << "[PHASE 2] INITIAL GEOMETRY GENERATION" << std::endl;
    std::cout << "------------------------------------------------------------------------" << std::endl;

    std::cout << "  [*] Setting up the initial simulation domain..." << std::endl;
    auto domain = viennaps::Domain<double, DIM>::New(
        params.get("gridDelta"),    
        params.get("xExtent"),      
        viennaps::BoundaryType::REFLECTIVE_BOUNDARY
    );

    std::cout << "  [*] Applying initial trench geometry..." << std::endl;
    viennaps::MakeTrench<double, DIM>(
        domain, 
        params.get("trenchWidth"),                
        params.get("trenchHeight"), 
        params.get("taperAngle"),
        params.get("maskHeight")
    ).apply();

    std::cout << "  [+] Initial geometry built.\n" << std::endl;



    // ========================================================================
    // [PHASE 3] THE OUTER LOOP: PROFILE EVOLUTION
    // ========================================================================
    int num_etch_steps = params.get<int>("numEtchSteps");         
    double etch_time_step = params.get("etchTimeStep");    
    double etch_yield_factor = params.get("etchYieldFactor"); 
    int num_charge_iterations = params.get("numChargeIterations");
    int num_particles = params.get<int>("numRaysPerPoint") * 10;
    double grid_delta = params.get("gridDelta");

    for (int etch_step = 0; etch_step < num_etch_steps; ++etch_step) {

        std::cout << "========================================================================" << std::endl;
        std::cout << " [PHASE 3] ETCH STEP " << (etch_step + 1) << " / " << num_etch_steps << std::endl;
        std::cout << "========================================================================" << std::endl;

        // --------------------------------------------------
        // [3.1] EXTRACT CURRENT GEOMETRY
        // --------------------------------------------------
        std::cout << "\n  [3.1] EXTRACTING VOLUMETRIC AND SURFACE DATA" << std::endl;
        std::cout << "  --------------------------------------------------" << std::endl;
        
        // Create a fresh temporary domain to prevent CellSet accumulation bugs
        auto extraction_domain = viennaps::Domain<double, DIM>::New();
        for (auto ls : domain->getLevelSets()) {
            extraction_domain->insertNextLevelSetAsMaterial(ls, viennaps::Material::Si, false);
        }


        // Re-generate the 2D grid and 1D mesh based on the NEW level-set positions
        extraction_domain->generateCellSet(
            params.get("maskHeight") + 10 * params.get("gridDelta"), 
            viennaps::Material::GAS, 
            true    
        );
        auto cellSet = extraction_domain->getCellSet();
        cellSet->buildNeighborhood();
        cellSet->addScalarData("ChargeDensity", 0.0);
        cellSet->addScalarData("ElectricPotential", 0.0);
        cellSet->addScalarData("E_x", 0.0);
        cellSet->addScalarData("E_y", 0.0);

        auto ExPtr = cellSet->getScalarData("E_x");
        auto EyPtr = cellSet->getScalarData("E_y");
        auto chargeDensityPtr = cellSet->getScalarData("ChargeDensity");
        int num_cells = cellSet->getElements().size();
        double cell_volume = grid_delta * grid_delta; 

        auto surfaceMesh = extraction_domain->getDiskMesh(); 
        std::string initial_trench_name = "01_trench_surface_step" + std::to_string(etch_step) + ".vtp";
        viennals::VTKWriter<double>(surfaceMesh, initial_trench_name).apply();
        std::cout << "  [*] Geometry extracted and saved to: '" << initial_trench_name << "'" << std::endl;

        // Build Embree Scene
        GeometryCollider<DIM> collider;
        collider.buildSurface(surfaceMesh, grid_delta);

        // Domain & Kill Bounds for Particle Tracing
        double current_min_y = params.get("maskHeight") + 10.0;
        for (const auto& node : surfaceMesh->nodes) {
            if (node[1] < current_min_y) {
                current_min_y = node[1];
            }
        }

        double x_extent = params.get("xExtent");
        double start_y = params.get("maskHeight") + 10.0;
        double y_bottom_kill_bound = current_min_y - 20.0; 
        double y_top_kill_bound = start_y; 
        double x_kill_bound = x_extent / 2.0;
        double theta_max = M_PI / 2.0; 

        std::mt19937 gen(std::random_device{}());
        std::uniform_real_distribution<double> dist_x(-(x_extent / 2.0), (x_extent / 2.0));
        auto normals = surfaceMesh->getCellData().getVectorData("Normals");

        std::vector<double> surfaceCharge(surfaceMesh->nodes.size(), 0.0);
        std::vector<double> ionFlux(surfaceMesh->nodes.size(), 0.0);
        std::vector<double> ionEnergy(surfaceMesh->nodes.size(), 0.0);
        std::vector<double> ionAngle(surfaceMesh->nodes.size(), 0.0);
        
        // Insert the array into the mesh once before the loop starts
        surfaceMesh->getPointData().insertNextScalarData(surfaceCharge, "SurfaceCharge");
        surfaceMesh->getPointData().insertNextScalarData(ionFlux, "IonFlux");
        surfaceMesh->getPointData().insertNextScalarData(ionEnergy, "IonEnergy");
        surfaceMesh->getPointData().insertNextScalarData(ionAngle, "IonAngle");

        // Pointers to the specific arrays inside the mesh memory
        auto exportedChargePtr = surfaceMesh->getPointData().getScalarData("SurfaceCharge");
        auto exportedFluxPtr = surfaceMesh->getPointData().getScalarData("IonFlux");
        auto exportedEnergyPtr = surfaceMesh->getPointData().getScalarData("IonEnergy");
        auto exportedAnglePtr = surfaceMesh->getPointData().getScalarData("IonAngle");



        // ========================================================================
        // [3.2] THE INNER LOOP: ELECTROSTATIC STEADY STATE
        // ========================================================================
        std::cout << "\n  [3.2] COMPUTING ELECTROSTATIC STEADY STATE" << std::endl;
        std::cout << "  --------------------------------------------------" << std::endl;

        for (int iter = 0; iter < num_charge_iterations; ++iter) {
            std::cout << "\n      >>> CHARGE ITERATION " << (iter + 1) << " / " << num_charge_iterations << " <<<" << std::endl;
            
            // --- Physics & Kinematics Setup ---
            double V_bias = params.get("V_bias");
            double d_sh = params.get("d_sh");
            double E_ext = V_bias / d_sh; 
            
            // Ion Parameters
            double mass_ion = params.get("mass");
            double q_ion = 1.0;
            double v_ion = std::sqrt(2.0 * q_ion * E_ext * d_sh / mass_ion);
            double n_theta_ion = 100.0; 
            double dt_ion = 0.01; 
            
            // Electron Parameters
            double mass_ratio = 1836.0; 
            double mass_elec = mass_ion / mass_ratio;
            double q_elec = -1.0;
            double v_elec = v_ion * std::sqrt(mass_ratio); 
            double n_theta_elec = 1.0;  
            double dt_elec = dt_ion / std::sqrt(mass_ratio); 

            // Macro-weights & Visualization Limits
            double macro_charge_weight_ion = params.get("ionWeight"); 
            double macro_charge_weight_elec = params.get("elecWeight");
            int particles_to_export = std::min(num_particles, 30); 
            
            std::vector<Particle> exported_ions; 
            std::vector<Particle> exported_elecs; 
            std::vector<int> ion_hits(surfaceMesh->nodes.size(), 0);
            std::vector<int> elec_hits(surfaceMesh->nodes.size(), 0);
            std::vector<double> ion_energy_sum(surfaceMesh->nodes.size(), 0);
            std::vector<double> elec_energy_sum(surfaceMesh->nodes.size(), 0);
            std::vector<double> ion_angle_sum(surfaceMesh->nodes.size(), 0);
            std::vector<double> elec_angle_sum(surfaceMesh->nodes.size(), 0);


            // --------------------------------------------------
            // ION TRACING
            // --------------------------------------------------
            std::cout << "  [*] Launching a batch of " << num_particles << " ions..." << std::endl;
            int hit_count_ion = 0;

            for (int i = 0; i < num_particles; ++i) {
                double start_x = dist_x(gen);
                double theta = samplePowerCosineAngle(gen, n_theta_ion, theta_max);
                double v_x = v_ion * std::sin(theta);
                double v_y = -v_ion * std::cos(theta); 

                Particle particle(mass_ion, q_ion, {start_x, start_y, 0.0}, {v_x, v_y, 0.0});
                bool active = true;

                while (active) {
                    std::array<double, 3> old_pos = particle.position;
                    double E_x = 0.0;
                    double E_y = 0.0;

                    int cell_idx = cellSet->getIndex(particle.position);
                    if (cell_idx >= 0) {
                        E_x = ExPtr->at(cell_idx);
                        E_y = EyPtr->at(cell_idx);
                    }

                    particle.advanceTimeStep(dt_ion, E_x, E_y);

                    double dx = particle.position[0] - old_pos[0];
                    double dy = particle.position[1] - old_pos[1];
                    double step_distance = sqrt(dx*dx + dy*dy);

                    // Only check for collisions if the particle actually moves, avoiding div0
                    if (step_distance > 1e-12) {
                        std::array<double, 3> direction = {dx/step_distance, dy/step_distance, 0.0};
                        int hit_id = collider.checkCollision(old_pos, direction, step_distance);
                        if (hit_id >= 0) {
                            ion_hits[hit_id]++;
                            ion_energy_sum[hit_id] += particle.energy;

                            // Calculate the angle of incidence relative to the surface normal
                            std::array<double, 3> normal = normals->at(hit_id);
                            double v_mag = std::sqrt(particle.velocity[0]*particle.velocity[0] + particle.velocity[1]*particle.velocity[1]);

                            // Dot product of velocity and normal. (Normal points OUT, velocity points IN, so we negate)
                            double dot_prod = -(particle.velocity[0] * normal[0] + particle.velocity[1] * normal[1]);
                            double cos_theta = std::max(-1.0, std::min(1.0, dot_prod / v_mag)); // Clamp to prevent NaN

                            ion_angle_sum[hit_id] += std::acos(cos_theta);
                            hit_count_ion++;
                            active = false; 
                        }
                    }
                    // Bounds checking (Kill particle if it flies out of the domain)
                    if (particle.position[1] < y_bottom_kill_bound || 
                        particle.position[1] > y_top_kill_bound || 
                        std::abs(particle.position[0]) > x_kill_bound) {
                        
                        active = false; 
                    }
                }
                if (i < particles_to_export) {
                    exported_ions.push_back(particle);
                }
            }


            // --------------------------------------------------
            // ELECTRON TRACING
            // --------------------------------------------------
            std::cout << "  [*] Launching a batch of " << num_particles << " electrons..." << std::endl;
            int hit_count_elec = 0;

            for (int i = 0; i < num_particles; ++i) {
                double start_x = dist_x(gen);
                double theta = samplePowerCosineAngle(gen, n_theta_elec, theta_max);
                double v_x = v_elec * std::sin(theta);
                double v_y = -v_elec * std::cos(theta); 

                Particle particle(mass_elec, q_elec, {start_x, start_y, 0.0}, {v_x, v_y, 0.0});
                bool active = true;

                while (active) {
                    std::array<double, 3> old_pos = particle.position;
                    double E_x = 0.0;
                    double E_y = 0.0;

                    int cell_idx = cellSet->getIndex(particle.position);
                    if (cell_idx >= 0) {
                        E_x = ExPtr->at(cell_idx);
                        E_y = EyPtr->at(cell_idx);
                    }

                    particle.advanceTimeStep(dt_elec, E_x, E_y);

                    double dx = particle.position[0] - old_pos[0];
                    double dy = particle.position[1] - old_pos[1];

                    double step_distance = sqrt(dx*dx + dy*dy);

                    // Only check for collisions if the particle actually moves, avoiding div0
                    if (step_distance > 1e-12) {
                        std::array<double, 3> direction = {dx/step_distance, dy/step_distance, 0.0};
                        int hit_id = collider.checkCollision(old_pos, direction, step_distance);
                        if (hit_id >= 0) {
                            elec_hits[hit_id]++;
                            elec_energy_sum[hit_id] += particle.energy;

                            // Calculate the angle of incidence relative to the surface normal
                            std::array<double, 3> normal = normals->at(hit_id);
                            double v_mag = std::sqrt(particle.velocity[0]*particle.velocity[0] + particle.velocity[1]*particle.velocity[1]);
                            
                            // Dot product of velocity and normal. (Normal points OUT, velocity points IN, so we negate)
                            double dot_prod = -(particle.velocity[0] * normal[0] + particle.velocity[1] * normal[1]);
                            double cos_theta = std::max(-1.0, std::min(1.0, dot_prod / v_mag)); // Clamp to prevent NaN

                            elec_angle_sum[hit_id] += std::acos(cos_theta);
                            hit_count_elec++;
                            active = false; 
                        }
                    }
                    // Bounds checking (Kill particle if it flies out of the domain)
                    if (particle.position[1] < y_bottom_kill_bound || 
                        particle.position[1] > y_top_kill_bound || 
                        std::abs(particle.position[0]) > x_kill_bound) {
                        
                        active = false; 
                    }
                }
                if (i < particles_to_export) {
                    exported_elecs.push_back(particle);
                }
            }

            std::cout << "              -> Recorded Hits: " << hit_count_ion << " Ions | " << hit_count_elec << " Electrons." << std::endl;


            // --------------------------------------------------
            // CHARGE CONTINUITY & AVERAGING
            // --------------------------------------------------
            std::cout << "          [*] Updating Surface Charge via Explicit Continuity..." << std::endl;
            // Pseudo-Time step (move into Parameters.txt?)
            double dt_c = params.get("dt_c"); 

            for (size_t i = 0; i < surfaceMesh->nodes.size(); ++i) {
                // Calculate Averages for Etch Model
                ionFlux[i] = static_cast<double>(ion_hits[i]) / num_particles;
                if (ion_hits[i] > 0) {
                    ionEnergy[i] = ion_energy_sum[i] / ion_hits[i];
                    ionAngle[i] = ion_angle_sum[i] / ion_hits[i];
                } else {
                    ionEnergy[i] = 0.0;
                    ionAngle[i] = 0.0;
                }
            }

            // Calculate Net Charge Accumulation
            for (size_t i = 0; i < surfaceCharge.size(); ++i) {
                // Convert raw hits to dimensionless normalized fluxes (hits / total launched)
                double gamma_i = ionFlux[i];
                double gamma_e = static_cast<double>(elec_hits[i]) / num_particles;

                // Calculate the net charge rate taking into account particle charge and macroscopic weighting
                double net_charge_rate = (q_ion * gamma_i * macro_charge_weight_ion) + 
                                        (q_elec * gamma_e * macro_charge_weight_elec);
                
                // Advance the explicit time step: Sigma_new = Sigma_old + dt * rate
                surfaceCharge[i] += dt_c * net_charge_rate;
            }

            // Bind data back to pointers for ParaView export and Advection reading
            *exportedChargePtr = surfaceCharge;
            *exportedFluxPtr = ionFlux;
            *exportedEnergyPtr = ionEnergy;
            *exportedAnglePtr = ionAngle;


            // --------------------------------------------------
            // PARAVIEW EXPORT
            // --------------------------------------------------
            std::string charge_filename = "02_charge_step" + std::to_string(etch_step) + "_iter" + std::to_string(iter) + ".vtp";
            viennals::VTKWriter<double>(surfaceMesh, charge_filename).apply();

            auto pathMeshIon = viennals::Mesh<double>::New();
            for (const auto& p : exported_ions) {
                unsigned int current_nodes = pathMeshIon->nodes.size();
                for (size_t j = 0; j < p.trajectory.size(); ++j) {
                    pathMeshIon->nodes.push_back(p.trajectory[j]);
                    if (j > 0) {
                        pathMeshIon->insertNextLine({static_cast<unsigned int>(current_nodes + j - 1), static_cast<unsigned int>(current_nodes + j)});
                    }
                }
            }
            std::string path_filename_ion = "03a_ion_paths_step" + std::to_string(etch_step) + "_iter" + std::to_string(iter) + ".vtp";
            viennals::VTKWriter<double>(pathMeshIon, path_filename_ion).apply();

            auto pathMeshElec = viennals::Mesh<double>::New();
            for (const auto& p : exported_elecs) {
                unsigned int current_nodes = pathMeshElec->nodes.size();
                for (size_t j = 0; j < p.trajectory.size(); ++j) {
                    pathMeshElec->nodes.push_back(p.trajectory[j]);
                    if (j > 0) {
                        pathMeshElec->insertNextLine({static_cast<unsigned int>(current_nodes + j - 1), static_cast<unsigned int>(current_nodes + j)});
                    }
                }
            }
            std::string path_filename_elec = "03b_elec_paths_step" + std::to_string(etch_step) + "_iter" + std::to_string(iter) + ".vtp";
            viennals::VTKWriter<double>(pathMeshElec, path_filename_elec).apply();


            // --------------------------------------------------
            // POISSON SOLVER
            // --------------------------------------------------
            std::cout << "          [*] Solving Poisson Equation for Electric Potential..." << std::endl;
            chargeDensityPtr->assign(num_cells, 0.0);

            for (size_t i = 0; i < surfaceCharge.size(); ++i) {
                if (std::abs(surfaceCharge[i]) > 1e-12) {
                    std::array<double, 3> node_coord = surfaceMesh->nodes[i];
                    std::array<double, 3> normal = normals->at(i);

                    node_coord[0] += normal[0] * (grid_delta * 0.5);
                    node_coord[1] += normal[1] * (grid_delta * 0.5);

                    int cell_idx = cellSet->getIndex(node_coord);
                    if (cell_idx >= 0) {
                        (*chargeDensityPtr)[cell_idx] += surfaceCharge[i] / cell_volume;
                    }
                }
            }
            PoissonSolver::solvePoisson(extraction_domain, grid_delta);

            std::string vol_filename = "04_volume_poisson_step" + std::to_string(etch_step) + "_iter" + std::to_string(iter) + ".vtu";
            extraction_domain->getCellSet()->writeVTU(vol_filename);
        }

        // ========================================================================
        // [3.3] ADVECTION (APPLY THE ETCH)
        // ========================================================================
        std::cout << "\n  [3.3] APPLYING SURFACE ETCH" << std::endl;
        std::cout << "  --------------------------------------------------" << std::endl;
        std::cout << "  [*] Steady State Reached. Initiating ViennaLS Advection..." << std::endl;
        
        // Pass the final, stabilized physics data into the Yield Model
        auto velField = viennals::SmartPointer<EtchVelocityField>(
            new EtchVelocityField(&(surfaceMesh->nodes), &ionFlux, &ionEnergy, &ionAngle, etch_yield_factor)
        );

        viennals::Advect<double, DIM> advection;

        // Move all solid materials in the domain
        for (auto ls : domain->getLevelSets()) {
            advection.insertNextLevelSet(ls);
        }
        advection.setVelocityField(velField);
        advection.setAdvectionTime(etch_time_step);
        
        advection.apply();
        std::cout << "  [+] Advection complete. Profile successfully evolved!\n" << std::endl;
    }

    std::cout << "========================================================================" << std::endl;
    std::cout << "                   SIMULATION COMPLETE SUCCESSFULLY                     " << std::endl;
    std::cout << "========================================================================\n" << std::endl;
    
    return 0;
}