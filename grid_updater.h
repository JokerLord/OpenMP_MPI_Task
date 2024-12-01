#pragma once
#include <cmath>
#include <utility>

static const double pi = atan(1) * 4;

class WaveEquationSolution {
public:
    const double Lx, Ly, Lz;
    const double aT;

    WaveEquationSolution(double Lx, double Ly, double Lz)
        : Lx(Lx), Ly(Ly), Lz(Lz),
          aT(pi * sqrt(4.0 / (Lx * Lx) + 4.0 / (Ly * Ly) + 1.0 / (Lz * Lz))) {}

    double calcValue(double x, double y, double z, double t) const {
        return sin((2 * pi * x / Lx) + 3 * pi) * sin((2 * pi * y / Ly) + 2 * pi) * sin(pi * z / Lz) * cos((aT * t) + pi);
    }
};


class GridUpdater {
    const WaveEquationSolution waveFunction;
    const double timeStep;

public:
    GridUpdater(double Lx, double Ly, double Lz, double timeStep) : waveFunction(Lx, Ly, Lz), timeStep(timeStep) {}

    double analyticalValue(const Grid &grid, size_t i, size_t j, size_t k, double t) const {
        double x = i * grid.stepX;
        double y = j * grid.stepY;
        double z = k * grid.stepZ;
        return waveFunction.calcValue(x, y, z, t);
    }

    void initializeGridT0(Grid &gridT0) const {
        #pragma omp parallel for collapse(3)
        for (size_t i = 0; i < gridT0.dimX; ++i) {
            for (size_t j = 0; j < gridT0.dimY; ++j) {
                for (size_t k = 0; k < gridT0.dimZ; ++k) {
                    gridT0.setValue(i, j, k, analyticalValue(gridT0, i, j, k, 0));
                }
            }
        }
    }

    double initializeGridT1(Grid &gridT0, Grid &gridT1) const {
        double maxError = 0.0;

        #pragma omp parallel for collapse(3) reduction(max:maxError)
        for (size_t i = 1; i < gridT1.dimX - 1; ++i) {
            for (size_t j = 1; j < gridT1.dimY - 1; ++j) {
                for (size_t k = 1; k < gridT1.dimZ -1; ++k) {
                    double uT0 = gridT0.getValue(i, j, k);
                    double uT1 = uT0 + (timeStep * timeStep / 2.0) * gridT0.calcLaplace(i, j, k);
                    gridT1.setValue(i, j, k, uT1);

                    double error = std::fabs(uT1 - analyticalValue(gridT1, i, j, k, timeStep));

                    maxError = std::max(maxError, error);
                }
            }
        }

        return maxError;
    }

    double updateGrid(Grid &prevGrid, Grid &curGrid, Grid &nextGrid, double t) const {
        double maxError = 0.0;

        #pragma omp parallel for collapse(3) reduction(max:maxError)
        for (size_t i = 0; i < nextGrid.dimX; ++i) {
            for (size_t j = 0; j < nextGrid.dimY; ++j) {
                for (size_t k = 1; k < nextGrid.dimZ - 1; ++k) {
                    double newValue = 2 * curGrid.getValue(i, j, k) - prevGrid.getValue(i, j, k) + timeStep * timeStep * curGrid.calcLaplace(i, j, k);
                    nextGrid.setValue(i, j, k, newValue);

                    double error = std::fabs(newValue - analyticalValue(nextGrid, i, j, k, t));

                    maxError = std::max(maxError, error);
                }
            }
        }

        return maxError;
    }

    void move(Grid &dst, Grid &src) {
        std::swap(dst, src);
    }
};
