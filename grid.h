#pragma once
#include <vector>

class Grid {
    double *values;

public:
    size_t dimX, dimY, dimZ;
    double stepX, stepY, stepZ;

    Grid(size_t pointNum, double Lx, double Ly, double Lz)
        : dimX(pointNum), dimY(pointNum), dimZ(pointNum),
          stepX(Lx / (pointNum - 1)), stepY(Ly / (pointNum - 1)), stepZ(Lz / (pointNum - 1)) {
            values = new double[dimX * dimY * dimZ];
          }

    size_t calcIndex(size_t i, size_t j, size_t k) const {
        return i * dimY * dimZ + j * dimZ + k;
    }

    double getValue(size_t x, size_t y, size_t z) const {
        return values[calcIndex(x, y, z)];
    }

    void setValue(size_t x, size_t y, size_t z, double value) {
        values[calcIndex(x, y, z)] = value;
    }

    double calcLaplace(size_t i, size_t j, size_t k) const {
        bool isBorderI = (i == 0 || i == dimX - 1);
        double dX = (getValue((isBorderI ? dimX - 1 : i) - 1, j, k)
                     - 2 * getValue(i, j, k)
                     + getValue((isBorderI ? 0 : i) + 1, j, k)) / (stepX * stepX);

        bool isBorderJ = (j == 0 || j == dimY - 1);
        double dY = (getValue(i, (isBorderJ ? dimY - 1 : j) - 1, k)
                     - 2 * getValue(i, j, k)
                     + getValue(i, (isBorderJ ? 0 : j) + 1, k)) / (stepY * stepY);

        double dZ = (getValue(i, j, k-1) - 2 * getValue(i, j, k) + getValue(i, j, k+1)) / (stepZ * stepZ);
        return dX + dY + dZ;
    }
};
