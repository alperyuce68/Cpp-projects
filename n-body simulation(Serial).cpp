#include <iostream> // Include the iostream library for input and output operations
#include <vector> // Include the vector library to use the std::vector container
#include <cmath> // Include the cmath library for mathematical functions such as std::sqrt
#include <chrono> // Include the chrono library for time-related functions
#include <thread> // Include the thread library for using sleep_for function

// Define a struct for a Particle
struct Particle {
    double x, y, z; // Position coordinates
    double vx, vy, vz; // Velocity components
    double mass; // Mass of the particle
};

// Function to compute forces between particles
void compute_forces(std::vector<Particle>& particles) {
    const double G = 6.67430e-11; // Gravitational constant
    for (size_t i = 0; i < particles.size(); ++i) { // Loop over each particle
        particles[i].vx = particles[i].vy = particles[i].vz = 0.0; // Reset velocities
        for (size_t j = 0; j < particles.size(); ++j) { // Loop over each particle to compute force
            if (i != j) { // Avoid self-interaction
                double dx = particles[j].x - particles[i].x; // Difference in x-coordinates
                double dy = particles[j].y - particles[i].y; // Difference in y-coordinates
                double dz = particles[j].z - particles[i].z; // Difference in z-coordinates
                double dist = std::sqrt(dx * dx + dy * dy + dz * dz); // Calculate distance between particles
                if (dist > 1e-9) { // Avoid division by zero
                    double F = (G * particles[i].mass * particles[j].mass) / (dist * dist); // Calculate force magnitude
                    particles[i].vx += F * dx / dist; // Update velocity component in x-direction
                    particles[i].vy += F * dy / dist; // Update velocity component in y-direction
                    particles[i].vz += F * dz / dist; // Update velocity component in z-direction
                }
            }
        }
    }
}

// Function to update positions of particles
void update_positions(std::vector<Particle>& particles, double dt) {
    for (auto& p : particles) { // Loop over each particle
        p.x += p.vx * dt; // Update x-coordinate
        p.y += p.vy * dt; // Update y-coordinate
        p.z += p.vz * dt; // Update z-coordinate
    }
}

int main() {
    const int num_particles = 1000; // Number of particles in the simulation
    const double dt = 0.01; // Time step for the simulation
    std::vector<Particle> particles(num_particles); // Create a vector of particles

    // Initialize particles with positions closer together
    double pos_range = 0.1; // Reduced range for closer initial positions
    for (auto& p : particles) { // Loop over each particle
        p.x = (rand() / (double)RAND_MAX) * pos_range; // Initialize x-coordinate
        p.y = (rand() / (double)RAND_MAX) * pos_range; // Initialize y-coordinate
        p.z = (rand() / (double)RAND_MAX) * pos_range; // Initialize z-coordinate
        p.mass = rand() / (double)RAND_MAX + 1.0; // Initialize mass
    }

    auto start = std::chrono::high_resolution_clock::now(); // Start timer

    // Simulation loop
    const int num_steps = 100; // Number of simulation steps
    for (int step = 0; step < num_steps; ++step) { // Loop for each simulation step
        compute_forces(particles); // Compute forces between particles
        update_positions(particles, dt); // Update positions based on velocities

        // Print positions of the first 10 particles (or adjust as needed)
        std::cout << "Step " << step << ":\n";
        for (int i = 0; i < 10 && i < num_particles; ++i) { // Loop over the first 10 particles
            std::cout << "Particle " << i << ": ("
                << particles[i].x << ", "
                << particles[i].y << ", "
                << particles[i].z << ")\n"; // Print particle position
        }
        std::cout << "\n";

        // Sleep to simulate real-time visualization
        std::this_thread::sleep_for(std::chrono::milliseconds(100)); // Adjust delay as needed
    }

    auto end = std::chrono::high_resolution_clock::now(); // End timer
    std::chrono::duration<double> elapsed = end - start; // Calculate elapsed time
    std::cout << "Serial execution time: " << elapsed.count() << " seconds\n"; // Print elapsed time

    return 0; // Return success
}
