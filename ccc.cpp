#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <conio.h>
#include <omp.h>
#include <chrono>  
#include <algorithm>
const float PI = 3.14159265359f;
const float G = 0.05;
const float smoothingRadius = 0.01f;
const float targetDensity = 0.5f;
const float pressureMultiplier = 2000.5f;
const float gravityMultiplier = -0.06f;
const float deltaTime = 0.03f;
const float dampingFactor = 3.0f;

struct Particle {
    float x, y;
    float vx, vy;
    float density, pressure;
    int cloud;
};

std::vector<Particle> particles;

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

void computeDensities() {
    #pragma omp parallel for
    for (size_t i = 0; i < particles.size(); ++i) {
        particles[i].density = 0;
        float xi = particles[i].x;
        float yi = particles[i].y;

        for (size_t j = 0; j < particles.size(); ++j) {
            float dx = particles[j].x - xi;
            float dy = particles[j].y - yi;
            float dst = std::sqrt(dx * dx + dy * dy);
            particles[i].density += smoothingKernel(dst, smoothingRadius);
        }

        particles[i].pressure = pressureMultiplier * (particles[i].density - targetDensity);
    }
}

void computeForces() {
    const float minDistance = 0.01f;       // Softening parameter to avoid infinite forces
    const float maxAcceleration = 20.0f;  // Cap to prevent unphysical acceleration

    #pragma omp parallel for
    for (size_t i = 0; i < particles.size(); ++i) {
        float xi = particles[i].x;
        float yi = particles[i].y;
        float ax = 0, ay = 0;

        for (size_t j = 0; j < particles.size(); ++j) {
            if (i == j) continue;

            float dx = particles[j].x - xi;
            float dy = particles[j].y - yi;
            float dst = std::sqrt(dx * dx + dy * dy);

            // Avoid zero division
            float effectiveDst = std::max(dst, minDistance);

            float dirX = dx / effectiveDst;
            float dirY = dy / effectiveDst;

            // Gravity with softening
            float softenedDst2 = effectiveDst * effectiveDst + 0.01f;  // Small softening
            ax -= gravityMultiplier * G * dirX / softenedDst2;
            ay -= gravityMultiplier * G * dirY / softenedDst2;

            // Pressure
            float sharedPressure = (particles[i].pressure / (particles[i].density * particles[i].density)) +
                                   (particles[j].pressure / (particles[j].density * particles[j].density));
            float smoothingDerivative = smoothingKernelDerivative(effectiveDst, smoothingRadius);
            ax -= dirX * sharedPressure * smoothingDerivative;
            ay -= dirY * sharedPressure * smoothingDerivative;

            // Damping
            if (effectiveDst < smoothingRadius) {
                float rvx = particles[i].vx - particles[j].vx;
                float rvy = particles[i].vy - particles[j].vy;
                float relativeVelocity = rvx * dirX + rvy * dirY;

                if (relativeVelocity < 0) {
                    float dampingForce = dampingFactor * relativeVelocity;
                    ax -= dampingForce * dirX;
                    ay -= dampingForce * dirY;
                }
            }
        }
        
        ax = std::clamp(ax, -maxAcceleration, maxAcceleration);
        ay = std::clamp(ay, -maxAcceleration, maxAcceleration);

        particles[i].vx += ax * deltaTime;
        particles[i].vy += ay * deltaTime;
    }
}


void updatePositions() {
    #pragma omp parallel for
    for (size_t i = 0; i < particles.size(); ++i) {
        particles[i].x += particles[i].vx * deltaTime;
        particles[i].y += particles[i].vy * deltaTime;
    }
}

void initializeParticles(int numLength = 25, int numWidth = 5) {
    float spacing = 20.0f;
    float armAngle = PI / 4.0f;
    float overlapThreshold = spacing * 0.8f;

    // V-shaped cloud
    for (int i = 0; i < numLength / 1.5; ++i) {
        for (int j = -numWidth; j <= numWidth; ++j) {
            float along = i * spacing;
            float offset = j * spacing;
            float randomOffset = (rand() % 100 - 50) * 0.1f;

            // Left arm before rotation
            float xl = -along * std::cos(armAngle) - (offset + randomOffset) * std::sin(armAngle);
            float yl = along * std::sin(armAngle) - (offset + randomOffset) * std::cos(armAngle);

            // Right arm before rotation
            float xr = along * std::cos(armAngle) + (offset + randomOffset) * std::sin(armAngle);
            float yr = along * std::sin(armAngle) - (offset + randomOffset) * std::cos(armAngle);

            // Distance between left and right arm (before rotation)
            float dx = xl - xr;
            float dy = yl - yr;
            float distance = std::sqrt(dx * dx + dy * dy);

            // Rotate both arms
            float xl_rot = yl;
            float yl_rot = -xl;
            float xr_rot = yr;
            float yr_rot = -xr;

            // Only add left-arm particle if it's not overlapping
            if (distance > overlapThreshold) {
                particles.push_back({xl_rot, yl_rot, -0.09f, -0.1f, 0.5f, 0.0f, 1});
            }

            // Always add right-arm particle
            particles.push_back({xr_rot, yr_rot, 0.2f, 0.1f, 0.5f, 0.0f, 1});
        }
    }

    // Rod-shaped cloud (unchanged)
    int rodLength = 50;
    for (int i = 0; i < 2*rodLength / 1.5; ++i) {
        for (int j = 0; j < 16; ++j) {
            float x = -100.0f + j * spacing;
            float y = 200.0f + i * spacing;
            float randomOffsetX = (rand() % 100 - 50) * 0.1f;
            float randomOffsetY = (rand() % 100 - 50) * 0.1f;

            x += randomOffsetX;
            y += randomOffsetY;

            float angle = -PI / 2.0f;
            float xr = x * std::cos(angle) - y * std::sin(angle);
            float yr = x * std::sin(angle) + y * std::cos(angle);

            particles.push_back({xr, yr, -0.2f, 0.02f, 0.5f, 0.0f, 2});
        }
    }
}


int main() {
    std::ofstream outFile("simulation_output.txt");
    if (!outFile.is_open()) {
        std::cerr << "Failed to open output file." << std::endl;
        return 1;
    }

    initializeParticles(50, 8);

    // Start timing
    auto start = std::chrono::high_resolution_clock::now();

    for (int step = 0; step < 600000; ++step) {
        computeDensities();
        computeForces();
        updatePositions();

        if (step % 10 == 0) {
            for (auto &p : particles) {
                outFile << p.x << " " << p.y << " " << p.cloud << "\n";
            }
        }
    }

    // End timing
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Simulation took " << elapsed.count() << " seconds.\n";

    outFile.close();
    getch();
    return 0;
    
}   