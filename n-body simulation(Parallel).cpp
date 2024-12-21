//Using the necessary libraries
#include <iostream> // Include the iostream library for input and output operations
#include <vector> // Include the vector library to use the std::vector container
#include <cmath> // Include the cmath library for mathematical functions such as std::sqrt
#include <chrono> // Include the chrono library for time-related functions
#include <thread> // Include the thread library for using sleep_for function
#include <omp.h>// Include the OpenMP library for parallel programming



// Particle structure definition
struct Particle {
    double x, y, z;   // Position of the particle in space
    double vx, vy, vz; // Velocity of the particle in space
    double mass;      // Mass of the particle
};

// Function to compute gravitational forces
void compute_forces(std::vector<Particle>& particles) {
    const double G = 6.67430e-11; // Gravitational constant
#pragma omp parallel for schedule(dynamic)// Parallelize the outer loop with dynamic scheduling

    for (size_t i = 0; i < particles.size(); ++i) {// Iterates over each particle to calculate forces
                                           
        particles[i].vx = particles[i].vy = particles[i].vz = 0.0; // Reset velocities

        for (size_t j = 0; j < particles.size(); ++j) { //Iterates over all particles 
                                                        //to compute interaction with particle i

            if (i != j) { // Avoid self-interaction
                double dx = particles[j].x - particles[i].x;
                double dy = particles[j].y - particles[i].y;
                double dz = particles[j].z - particles[i].z;
                double dist = std::sqrt(dx * dx + dy * dy + dz * dz);// Calculate distance
                if (dist > 1e-9) { // Avoid division by zero
                    double F = (G * particles[i].mass * particles[j].mass) / (dist * dist);// Compute gravitational force
                    particles[i].vx += F * dx / dist; // Update velocity in x direction
                    particles[i].vy += F * dy / dist; // Update velocity in y direction
                    particles[i].vz += F * dz / dist; // Update velocity in z direction
                }
            }
        }
    }
}


// Function to update particle positions based on velocities
void update_positions(std::vector<Particle>& particles, double dt) {
#pragma omp parallel for // Parallelize the loop
    for (auto& p : particles) { //Iterates over each particle to update its position
        p.x += p.vx * dt; // Update x position
        p.y += p.vy * dt; // Update y position
        p.z += p.vz * dt; // Update z position
    }
}

int main() { //main using for executing the our code or fuctiong
    const int num_particles = 1000; // Number of particles
    const double dt = 0.01; // Time step for the simulation
    std::vector<Particle> particles(num_particles); // Vector of particles

    // Initialize particles with positions closer together
    double pos_range = 0.01; // Reduced range for closer initial positions
    for (auto& p : particles) {//Initializes each particle with random positions and mass
        p.x = (rand() / (double)RAND_MAX) * pos_range; // Random x position within range
        p.y = (rand() / (double)RAND_MAX) * pos_range; // Random y position within range
        p.z = (rand() / (double)RAND_MAX) * pos_range; // Random z position within range
        p.mass = rand() / (double)RAND_MAX + 1.0; // Random mass between 1.0 and 2.0
    }

    // Start time measurement
    auto start = std::chrono::high_resolution_clock::now();


    // Simulation loop
    const int num_steps = 100; // Number of simulation steps
    for (int step = 0; step < num_steps; ++step) { //Simulation loop, iterates over each simulation step
        compute_forces(particles); // Compute forces between particles
        update_positions(particles, dt); // Update positions based on velocities


        // Print positions of the first 10 particles (or adjust as needed)
        std::cout << "Step " << step << ":\n";
        for (int i = 0; i < 10 && i < num_particles; ++i) {//Prints the positions of the first 10 particles
            std::cout << "Particle " << i << ": ("
                << particles[i].x << ", "
                << particles[i].y << ", "
                << particles[i].z << ")\n";// Output positions of first 10 particles
        }
        std::cout << "\n";

        // Sleep to simulate real-time visualization
        std::this_thread::sleep_for(std::chrono::milliseconds(100)); // Adjust delay as needed
    }
    
    auto end = std::chrono::high_resolution_clock::now();// End time measurement
    std::chrono::duration<double> elapsed = end - start;// Calculate elapsed time
    std::cout << "Parallel execution time: " << elapsed.count() << " seconds\n";// Output execution time

    return 0;// Indicate successful execution
}
