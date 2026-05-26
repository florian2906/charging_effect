#pragma once
#include <array>
#include <cmath>

class Ion {
    public:
    double mass;
    double charge;
    std::array<double, 3> position;
    std::array<double, 3> velocity;
    std::vector<std::array<double, 3>> trajectory; 

    Ion(double m, double q, const std::array<double, 3>& pos, const std::array<double, 3>& vel)
        : mass(m), charge(q), position(pos), velocity(vel) {}

    void advanceTimeStep(double dt, double E_x, double E_y) {
        velocity[0] += charge / mass * E_x * dt;
        velocity[1] += charge / mass * E_y * dt;

        position[0] += velocity[0] * dt;
        position[1] += velocity[1] * dt;

        trajectory.push_back(position);
    }
};