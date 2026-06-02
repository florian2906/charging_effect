#pragma once
#include <psDomain.hpp>
#include <vector>
#include <cmath>

namespace PoissonSolver {

    inline void solvePoisson(viennaps::SmartPointer<viennaps::Domain<double, 2>> domain, double grid_delta) {
        auto cellSet = domain->getCellSet();
        auto materialIDs = cellSet->getScalarData("Material");
        
        // The total number of cells in the simulation
        int num_cells = materialIDs->size(); 

        auto chargeDensityPtr = cellSet->getScalarData("ChargeDensity");
        if (!chargeDensityPtr) {
            std::cerr << "\n[!] FATAL: ChargeDensity array not found!" << std::endl;
            return;
        }

        auto potentialPtr = cellSet->getScalarData("ElectricPotential");
        if (!potentialPtr) {
            potentialPtr = cellSet->addScalarData("ElectricPotential", 0.0);
            potentialPtr->assign(num_cells, 0.0);
        }

        auto ExPtr = cellSet->getScalarData("E_x");
        if (!ExPtr) {
            ExPtr = cellSet->addScalarData("E_x", 0.0);
            ExPtr->assign(num_cells, 0.0);
        }

        auto EyPtr = cellSet->getScalarData("E_y");
        if (!EyPtr) {
            EyPtr = cellSet->addScalarData("E_y", 0.0);
            EyPtr->assign(num_cells, 0.0);
        }

        double dx2 = grid_delta * grid_delta;
        double eps_0 = 1.0;

        int numSteps = 1000;
        
        std::cout << "  [*] Starting Jacobi iterations for Poisson solver with " << num_cells << " cells." << std::endl;
        // ==========================================
        // 1. JACOBI ITERATION: Solve for Potential (Phi)
        // ==========================================
        double tolerance = 1e-6;

        for (int step = 0; step < numSteps; ++step) {
            std::vector<double> next_potential(*potentialPtr);
            double max_change = 0.0;

            std::cout << "      -> Iterating step " << step << " / " << numSteps << "\r" << std::flush;

            for (int i = 0; i < num_cells; ++i) {
                // only solve the field inside the GAS
                if (viennaps::MaterialMap::mapToMaterial(materialIDs->at(i)) != viennaps::Material::GAS) {
                    continue;
                }

                int numNeighbors = 0;
                double neighbor_sum = 0.0;
                auto cellNeighbors = cellSet->getNeighbors(i);

                for (int n : cellNeighbors) {
                    if (n != -1) { // valid neighbor
                        neighbor_sum += potentialPtr->at(n);
                        numNeighbors++;
                    }
                }

                // Central Difference Poisson Equation: grad^2(Phi) = -rho / eps_0
                if (numNeighbors > 0) {
                    double rho = chargeDensityPtr->at(i);
                    next_potential[i] = (neighbor_sum + (rho * dx2 / eps_0)) / static_cast<double>(numNeighbors);

                    double change = std::abs(next_potential[i] - potentialPtr->at(i));
                    if (change > max_change) {
                        max_change = change;
                    }
                }
            }

            *potentialPtr = next_potential; // Update the grid for the next sweep

            // Convergence Check: If the largest change is smaller than our tolerance, stop iterating
            if (max_change < tolerance) {
                std::cout << "\n      -> Poisson solver converged early at step " << step 
                          << " (Max change: " << max_change << ")" << std::endl;
                break; // Exit
            }        
        }

        std::cout << "\n  [*] Calculating Electric Field Gradient..." << std::endl;

        // ==========================================
        // 2. GRADIENT CALCULATION: E = -grad(Phi)
        // ==========================================
        for (int i = 0; i < num_cells; i++) {
            if (viennaps::MaterialMap::mapToMaterial(materialIDs->at(i)) != viennaps::Material::GAS) continue;

            auto cellNeighbors = cellSet->getNeighbors(i);
            auto leftN = cellNeighbors[0];
            auto rightN = cellNeighbors[1];
            auto downN = cellNeighbors[2];
            auto upN = cellNeighbors[3];

            // E_x = -(Phi_right - Phi_left) / 2dx
            if (leftN == -1 && rightN != -1) ExPtr->at(i) = -(potentialPtr->at(rightN) - potentialPtr->at(i)) / grid_delta;
            else if (rightN == -1 && leftN != -1) ExPtr->at(i) = -(potentialPtr->at(i) - potentialPtr->at(leftN)) / grid_delta;
            else if (leftN != -1 && rightN != -1) ExPtr->at(i) = -(potentialPtr->at(rightN) - potentialPtr->at(leftN)) / (2.0 * grid_delta);

            // E_y = -(Phi_up - Phi_down) / 2dx
            if (downN == -1 && upN != -1) EyPtr->at(i) = -(potentialPtr->at(upN) - potentialPtr->at(i)) / grid_delta;
            else if (upN == -1 && downN != -1) EyPtr->at(i) = -(potentialPtr->at(i) - potentialPtr->at(downN)) / grid_delta;
            else if (downN != -1 && upN != -1) EyPtr->at(i) = -(potentialPtr->at(upN) - potentialPtr->at(downN)) / (2.0 * grid_delta);
        }
    }

} // namespace PoissonSolver