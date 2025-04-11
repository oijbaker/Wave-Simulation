// import necessary libraries
#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <random>

// parameter set-up
const int x_dim = 70;
const int y_dim = 70;
const double c = 0.1;

// Mersenne twister random number generator
std::mt19937 rng(std::random_device{}());
std::uniform_real_distribution<double> unif(0.0, 1.0);

const std::string filename = "wave_output.txt";

// Set up excluded region as a circular central object
const int object_radius = 10; // radius of the circular object
const int object_center_x = x_dim / 2; // x-coordinate of the center of the object
const int object_center_y = y_dim / 2; // y-coordinate of the center of the object


bool check_collision(int x, int y) {
    // Check if the point (x, y) is inside the circular object
    return (std::pow(x - object_center_x, 2) + std::pow(y - object_center_y, 2)) <= std::pow(object_radius, 2);
}

// initial state set-up
void single_drop_init(std::vector<std::vector<double>>& u) {
    for (int i = 0; i < x_dim; ++i) {
        for (int j = 0; j < y_dim; ++j) {
            u[i][j] = 0.0;
        }
    }
    
    double drop_x = unif(rng);
    double drop_y = unif(rng);

    while (check_collision(drop_x * x_dim, drop_y * y_dim)) {
        drop_x = unif(rng);
        drop_y = unif(rng);
    }

    int drop_i = static_cast<int>(drop_x * x_dim);
    int drop_j = static_cast<int>(drop_y * y_dim);
    u[drop_i][drop_j] = 1.0; // initial drop in a random position
}

// update function
void update_wave(std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& u_prev, bool absorb) {
    std::vector<std::vector<double>> u_next(x_dim, std::vector<double>(y_dim, 0.0));
    for (int i = 1; i < x_dim - 1; ++i) {
        for (int j = 1; j < y_dim - 1; ++j) {
            u_next[i][j] = 2 * u[i][j] - u_prev[i][j] + c * c * (u[i + 1][j] + u[i - 1][j] + u[i][j + 1] + u[i][j - 1] - 4 * u[i][j]);
        }
    }
    u_prev = u;

    // Dirichlet boundary conditions
    if (absorb) {
        for (int i = 0; i < x_dim; ++i) {
            u_next[i][0] = u[i][1] - (c-1)/(c+1) * (u[i][0] - u_next[i][1]); // top boundary
            u_next[i][y_dim - 1] = u[i][y_dim - 2] - (c-1)/(c+1) * (u[i][y_dim - 1] - u_next[i][y_dim - 2]); // bottom boundary
        } 
        for (int j = 0; j < y_dim; ++j) {
            u_next[0][j] = u[1][j] - (c-1)/(c+1) * (u[0][j] - u_next[1][j]); // left boundary
            u_next[x_dim - 1][j] = u[x_dim - 2][j] - (c-1)/(c+1) * (u[x_dim - 1][j] - u_next[x_dim - 2][j]); // right boundary
        }
    } else {
        for (int i = 0; i < x_dim; ++i) {
            u_next[i][0] = 0.0; // top boundary
            u_next[i][y_dim - 1] = 0.0; // bottom boundary
        }
        for (int j = 0; j < y_dim; ++j) {
            u_next[0][j] = 0.0; // left boundary
            u_next[x_dim - 1][j] = 0.0; // right boundary
        }

        // check for collisions with the circular object
        for (int i = 0; i < x_dim; ++i) {
            for (int j = 0; j < y_dim; ++j) {
                if (check_collision(i, j)) {
                    u_next[i][j] = 0.0; // set the wave amplitude to zero inside the object
                }
            }
        }
    }
    u = u_next;
} 

void run(int steps, bool absorb) {
    std::vector<std::vector<double>> u(x_dim, std::vector<double>(y_dim, 0.0));
    std::vector<std::vector<double>> u_prev(x_dim, std::vector<double>(y_dim, 0.0));
    single_drop_init(u);

    std::cout << "Running simulation with " << (absorb ? "absorbing" : "non-absorbing") << " boundary conditions." << std::endl;

    // store all the states in a list
    std::vector<std::vector<std::vector<double>>> states(steps, std::vector<std::vector<double>>(x_dim, std::vector<double>(y_dim, 0.0)));
    for (int t = 0; t < steps; ++t) {
        states[t] = u;
        update_wave(u, u_prev, absorb);
        if (t % 10 == 0) { // print every 10 steps
            std::cout << (double) t/steps*100 << "% Complete" << std::endl;
        }
    }
    std::cout << "100% Complete" << std::endl;
    std::cout << "Writing output to file..." << std::endl;
    
    std::ofstream outfile(filename);
    if (!outfile) {
        std::cerr << "Error opening file for writing." << std::endl;
        return;
    }
    for (int t = 0; t < steps; ++t) {
        for (int i = 0; i < x_dim; ++i) {
            for (int j = 0; j < y_dim; ++j) {
                outfile << states[t][i][j] << " ";
            }
        }
        outfile << "\n";
    }
    outfile.close();

    std::cout << "Output written to " << filename << std::endl;
    std::cout << "Simulation finished." << std::endl;
}


int main() {
    int steps = 1000; // number of time steps to simulate
    bool absorb = false; // use absorbing boundary conditions
    run(steps, absorb);
    return 0;
}
