#include <cmath>
#include <fstream>
#include <iostream>
#include <cstdlib>  

const float PI = 3.14159265359f;
const float smoothingRadius = 1.0f;
const float targetDensity = 2.0f;
const float pressureMultiplier = 2.5f;
const float gravity = 90.8f;
const float deltaTime = 0.006f;
const int numParticles = 600;
const float dampingFactor = 0.85f;  
const float boxSize = 10.0f;  

float positionsX[numParticles];
float positionsY[numParticles];
float velocitiesX[numParticles] = {0};
float velocitiesY[numParticles] = {0};
float densities[numParticles] = {0};

float smoothingKernel(float dst, float radius) {
    if (dst >= radius) return 0;
    float volume = (PI * std::pow(radius, 4)) / 6.0f;
    return (radius - dst) * (radius - dst) / volume;
}

float smoothingKernelDerivative(float dst, float radius) {
    if (dst >= radius) return 0;
    float volume = (PI * std::pow(radius, 4)) / 6.0f;
    return -2.0f * (radius - dst) / volume;
}

float calculateDensity(float sampleX, float sampleY) {
    float density = 0;
    const float mass = 1.0f;
    for (int i = 0; i < numParticles; ++i) {
        float dx = positionsX[i] - sampleX;
        float dy = positionsY[i] - sampleY;
        float dst = std::sqrt(dx * dx + dy * dy);
        density += mass * smoothingKernel(dst, smoothingRadius);
    }
    return density;
}

float convertDensityToPressure(float density) {
    float densityError = density - targetDensity;
    return densityError * pressureMultiplier;
}

void applyBoundaryConditions() {
    for (int i = 0; i < numParticles; ++i) {
        // Reflective boundary conditions 
        if (positionsX[i] <= 0.0f && positionsY[i] != 0.0f){
            continue;
      }
        //if(positionsX[i] >= 2.5f && positionsX[i] <= 4.0f && positionsY[i] <= 0.0f ) {
          //  continue;
        //}
        if (positionsX[i] <= 0.0f) {           //
            positionsX[i] = 0.0f;              //
            velocitiesX[i] = -velocitiesX[i];   // Reflect velocity
        }                                      //
        
        if (positionsX[i] >= boxSize) {        // 
            positionsX[i] = boxSize;           //
            velocitiesX[i] = -velocitiesX[i] ;  //
        }                                      //
                                               //
        if (positionsY[i] <= 0.0f) {           //
            positionsY[i] = 0.0f;              //
            velocitiesY[i] = -velocitiesY[i]  ;//
        }                                      //
        if (positionsY[i] >= boxSize) {        //
            positionsY[i] = boxSize;           //
            velocitiesY[i] = -velocitiesY[i]  ;//
        }
    }
}

void simulationStep() {
    // Update densities
    for (int i = 0; i < numParticles; ++i) {
        densities[i] = calculateDensity(positionsX[i], positionsY[i]);
    }

    // Update velocities
    for (int i = 0; i < numParticles; ++i) {
        float pressureForceX = 0;
        float pressureForceY = 0;
        for (int j = 0; j < numParticles; ++j) {
            if (i == j) continue;
            float dx = positionsX[j] - positionsX[i];
            float dy = positionsY[j] - positionsY[i];
            float dst = std::sqrt(dx * dx + dy * dy);
            float dirX = dst > 0 ? dx / dst : 0;
            float dirY = dst > 0 ? dy / dst : 0;

            // Use the derivative of the smoothing kernel to calculate the pressure force
            float smoothingDerivative = smoothingKernelDerivative(dst, smoothingRadius);

            float densityB = densities[j];
            float sharedPressure = convertDensityToPressure(densities[i]) + convertDensityToPressure(densityB);

            pressureForceX += dirX * sharedPressure * smoothingDerivative;
            pressureForceY += dirY * sharedPressure * smoothingDerivative;
        }

        // Add gravity 
        velocitiesX[i] += pressureForceX * deltaTime;
        velocitiesY[i] += pressureForceY * deltaTime - gravity * deltaTime;

        // Apply damping
        velocitiesX[i] *= dampingFactor;
        velocitiesY[i] *= dampingFactor;
    }

    // Update positions
    for (int i = 0; i < numParticles; ++i) {
        positionsX[i] += velocitiesX[i] * deltaTime;
        positionsY[i] += velocitiesY[i] * deltaTime;
    }

    // Apply boundary conditions 
    applyBoundaryConditions();
}

void initializeParticles() {
    float range = 10.0f;  // Define the maximum coordinate range
    for (int i = 0; i < numParticles; ++i) {
        // Randomly place particles within the defined range
        positionsX[i] = static_cast<float>(rand()) / RAND_MAX * range/2;
        positionsY[i] = static_cast<float>(rand()) / RAND_MAX * range*2;
    }
}

int main() {
    std::ofstream outFile("simulation_output.txt");

    if (!outFile.is_open()) {
        std::cerr << "Failed to open output file." << std::endl;
        return 1;
    }

    initializeParticles();

    for (int step = 0; step < 3000; ++step) {
        simulationStep();
        if (step % 10 == 0) { // Log every 10 steps
            //outFile << "step " << step << ":\n";
            for (int i = 0; i < numParticles; ++i) {
                outFile << positionsX[i] << " " << positionsY[i] << "\n";
            }
            // outFile << "\n";
        }
    }

    outFile.close();
    return 0;
}