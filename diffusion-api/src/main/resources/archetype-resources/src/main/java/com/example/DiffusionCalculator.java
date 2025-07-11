package com.example;

public class DiffusionCalculator {

    public double calculateDiffusionCoefficient(double temperature, double pressure, double viscosity) {
        // This is a placeholder for your actual calculation logic
        // Replace this with the correct formula and parameters
        if (viscosity <= 0) {
            throw new IllegalArgumentException("Viscosity must be positive.");
        }
        return (1.86e-3 * Math.pow(temperature, 1.5)) / (pressure * viscosity);
    }
}