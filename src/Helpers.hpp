// Helper function to sample an angle from a power-cosine distribution P(theta) ~ cos^n(theta)

double samplePowerCosineAngle(std::mt19937& gen, double n_theta, double theta_max = M_PI / 2.0) {
    std::uniform_real_distribution<double> dist_theta(-theta_max, theta_max);
    std::uniform_real_distribution<double> dist_u(0.0, 1.0);
     // High n_theta value = tightly collimated vertical beam
     // Ions: n_theta = 100
     // Electrons: n_theta = 1
    while (true) {
        double theta = dist_theta(gen);
        double u = dist_u(gen);
        // Rejection condition
        if (u <= std::pow(std::cos(theta), n_theta)) {
            return theta;
        }
    }
}