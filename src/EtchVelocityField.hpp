#pragma once
#include <lsVelocityField.hpp>
#include <vector>
#include <array>
#include <cmath>
#include <limits>

class EtchVelocityField : public viennals::VelocityField<double> {
private:
    const std::vector<std::array<double, 3>>* surface_nodes;
    const std::vector<double>* flux;
    const std::vector<double>* energy;
    const std::vector<double>* angle;
    double etch_yield_factor;
    double mask_selectivity;

public:
    // Constructor takes pointers to our extracted mesh data so it can read the physics
    EtchVelocityField(const std::vector<std::array<double, 3>>* nodes, 
                      const std::vector<double>* f, 
                      const std::vector<double>* e, 
                      const std::vector<double>* a, 
                      double factor = 1.0,
                      double selectivity = 0.0)
        : surface_nodes(nodes), flux(f), energy(e), angle(a), etch_yield_factor(factor), mask_selectivity(selectivity) {}

    // This function is automatically called by ViennaLS for every grid point near the surface
    double getScalarVelocity(const std::array<double, 3>& coordinate, int material,
                             const std::array<double, 3>& normal, unsigned long pointID) override {
        
        // SPATIAL LOOKUP: Find the closest node on our Embree surface mesh 
        // to the current Level-Set grid coordinate being evaluated by ViennaLS.
        double min_dist = std::numeric_limits<double>::max();
        size_t closest_idx = 0;
        
        for (size_t i = 0; i < surface_nodes->size(); ++i) {
            double dx = coordinate[0] - (*surface_nodes)[i][0];
            double dy = coordinate[1] - (*surface_nodes)[i][1];
            double dist = dx*dx + dy*dy; // 2D distance squared
            
            if (dist < min_dist) {
                min_dist = dist;
                closest_idx = i;
            }
        }

        // READ PHYSICS: Get the steady-state plasma data for this specific node
        double f = (*flux)[closest_idx];
        double e = (*energy)[closest_idx];
        double a = (*angle)[closest_idx];

        // If no ions hit this node, it doesn't etch
        if (f <= 0.0) return 0.0;

        // MATERIAL SELECTIVITY
        // In ViennaPS, material 0 is the Substrate, material 1 is the Mask.
        double current_material_factor = (material == 1) ? mask_selectivity : 1.0;

        // PHYSICAL YIELD MODEL: 
        // Etch rate scales with the square root of energy and cosine of incident angle.
        double yield = std::sqrt(e) * std::cos(a);
        
        // Final surface velocity
        double etch_rate = f * yield * etch_yield_factor * current_material_factor;
        
        // Return a NEGATIVE velocity to tell the level-set to move "inward" (remove material)
        return -etch_rate;
    }

    // Required by ViennaLS interface, but we only use scalar (normal) velocity for etching
    std::array<double, 3> getVectorVelocity(const std::array<double, 3>& coordinate, int material,
                                            const std::array<double, 3>& normal, unsigned long pointID) override {
        return {0.0, 0.0, 0.0};
    }
};