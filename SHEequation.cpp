#include <iostream>  // Standard input/output stream
#include <cmath>     // Mathematical functions and constants
#include <vector>    // Dynamic array container
#include <iomanip>   // Input/output manipulators

#define M_PI 3.14159265358979323846264  // Definition of mathematical constant pi

// Function to define the SHE equations for a 2-level inverter
std::vector<double> she_equations(double x, double modulation_index, int order) {
    std::vector<double> equations(order, 0.0);

    // Calculate SHE equations for each harmonic
    for (int i = 1; i <= order; ++i) {
        equations[i - 1] = modulation_index * sin(i * M_PI * x) - x;
    }

    return equations;
}

// Function to solve SHE equations using the Newton-Raphson method
double newton(std::vector<double> (*equations)(double, double, int), double guess, double modulation_index, int order) {
    double x = guess;

    // Setting tolerance for convergence
    const double tolerance = 1e-6;

    // Maximum number of iterations
    const int max_iterations = 100;

    // Newton-Raphson iteration
    for (int iteration = 0; iteration < max_iterations; ++iteration) {
        std::vector<double> f = equations(x, modulation_index, order);

        // Checking for convergence
        if (std::abs(f[0]) < tolerance) {
            break;
        }

        // Updating x using Newton's method
        x -= f[0] / order;
    }

    return x;
}

// Main function
int main() {
    double modulation_index;
    int order;

    // Getting user input
    std::cout << "Enter modulation index: ";
    std::cin >> modulation_index;

    std::cout << "Enter the order of SHE equations: ";
    std::cin >> order;

    // Initial guess for Newton-Raphson
    double initial_guess = 0.5;

    // Solving SHE equations
    double solution = newton(she_equations, initial_guess, modulation_index, order);

    // Printing solution with 4 decimal places
    std::cout << "\nSHE Solution: x = " << std::fixed << std::setprecision(4) << solution << std::endl;

    return 0;
}
