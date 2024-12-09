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

    void initializeGridT0(Grid &gridT0, size_t *startIndices) const {
        #pragma omp parallel for collapse(3)
        for (size_t i = 1; i < gridT0.dimX - 1; ++i) {
            for (size_t j = 1; j < gridT0.dimY - 1; ++j) {
                for (size_t k = 1; k < gridT0.dimZ - 1; ++k) {
                    double analVal = analyticalValue(gridT0,
                                                     i + startIndices[0] - 1,
                                                     j + startIndices[1] - 1,
                                                     k + startIndices[2] - 1,
                                                     0
                                     );
                    gridT0.setValue(i, j, k, analVal);
                }
            }
        }
    }

    double initializeGridT1(Grid &gridT0, Grid &gridT1, size_t *startIndices) const {
        double maxError = 0.0;

        #pragma omp parallel for collapse(3) reduction(max:maxError)
        for (size_t i = 1; i < gridT1.dimX - 1; ++i) {
            for (size_t j = 1; j < gridT1.dimY - 1; ++j) {
                for (size_t k = 2; k < gridT1.dimZ - 2; ++k) {
                    double uT0 = gridT0.getValue(i, j, k);
                    double uT1 = uT0 + (timeStep * timeStep / 2.0) * gridT0.calcLaplace(i, j, k);
                    gridT1.setValue(i, j, k, uT1);

                    double analValue = analyticalValue(gridT1,
                                                        i + startIndices[0] - 1,
                                                        j + startIndices[1] - 1,
                                                        k + startIndices[2] - 1,
                                                        timeStep);

                    double error = std::fabs(uT1 - analValue);

                    maxError = std::max(maxError, error);
                }
            }
        }

        return maxError;
    }

    double updateGrid(Grid &prevGrid, Grid &curGrid, Grid &nextGrid, double t, size_t *startIndices) const {
        double maxError = 0.0;

        #pragma omp parallel for collapse(3) reduction(max:maxError)
        for (size_t i = 1; i < nextGrid.dimX - 1; ++i) {
            for (size_t j = 1; j < nextGrid.dimY - 1; ++j) {
                for (size_t k = 2; k < nextGrid.dimZ - 2; ++k) {
                    double newValue = 2 * curGrid.getValue(i, j, k) - prevGrid.getValue(i, j, k) + timeStep * timeStep * curGrid.calcLaplace(i, j, k);
                    nextGrid.setValue(i, j, k, newValue);

                    double analValue = analyticalValue(nextGrid,
                                                       i + startIndices[0] - 1,
                                                       j + startIndices[1] - 1,
                                                       k + startIndices[2] - 1,
                                                       t);

                    double error = std::fabs(newValue - analValue);

                    maxError = std::max(maxError, error);
                }
            }
        }

        return maxError;
    }

    void move(Grid &dst, Grid &src) {
        std::swap(dst.values, src.values);
    }
};
