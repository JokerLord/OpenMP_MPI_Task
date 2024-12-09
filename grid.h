#pragma once
#include <vector>

class Grid {
public:
    size_t globalPointNum;
    size_t dimX, dimY, dimZ;
    double stepX, stepY, stepZ;

    size_t boundaryDimX;
    size_t boundaryDimY;
    size_t boundaryDimZ;

    size_t *startIndices;

    double *values, *upperPart1, *lowerPart1, *frontPart1, *backPart1, *leftPart1, *rightPart1;
    double *upperPart2, *lowerPart2, *frontPart2, *backPart2, *leftPart2, *rightPart2;

    Grid(size_t pointNum, double Lx, double Ly, double Lz, size_t sizeX, size_t sizeY, size_t sizeZ, size_t *blockStartIndices) {
        globalPointNum = pointNum;

        dimX = sizeX;
        dimY = sizeY;
        dimZ = sizeZ;
        stepX = Lx / ((double) pointNum - 1);
        stepY = Ly / ((double) pointNum - 1);
        stepZ = Lz / ((double) pointNum - 1);
        values = new double[dimX * dimY * dimZ];

        upperPart1 = new double[dimX * dimZ];
        lowerPart1 = new double[dimX * dimZ];
        frontPart1 = new double[dimX * dimY];
        backPart1 = new double[dimX * dimY];
        leftPart1 = new double[dimY * dimZ];
        rightPart1 = new double[dimY * dimZ];

        upperPart2 = new double[dimX * dimZ];
        lowerPart2 = new double[dimX * dimZ];
        frontPart2 = new double[dimX * dimY];
        backPart2 = new double[dimX * dimY];
        leftPart2 = new double[dimY * dimZ];
        rightPart2 = new double[dimY * dimZ];

        startIndices = blockStartIndices;
    }

    size_t calcIndex(size_t i, size_t j, size_t k) const {
        return i * dimY * dimZ + j * dimZ + k;
    }

    size_t calcBoundaryIndex(size_t i, size_t j, size_t size) const {
        return i * size + j;
    }

    double getValue(size_t x, size_t y, size_t z) const {
        return values[calcIndex(x, y, z)];
    }

    double* getPtr(size_t x, size_t y, size_t z) {
        return &values[calcIndex(x, y, z)];
    }

    void setValue(size_t x, size_t y, size_t z, double value) {
        values[calcIndex(x, y, z)] = value;
    }

    double calcLaplace(size_t i, size_t j, size_t k) const {
        double dX = (getValue(i - 1, j, k) - 2 * getValue(i, j, k) + getValue(i + 1, j, k)) / (stepX * stepX);
        double dY = (getValue(i, j - 1, k) - 2 * getValue(i, j, k) + getValue(i, j + 1, k)) / (stepY * stepY);
        double dZ = (getValue(i, j, k - 1) - 2 * getValue(i, j, k) + getValue(i, j, k + 1)) / (stepZ * stepZ);
        return dX + dY + dZ;
    }

    void fillBoundary() {
        if (startIndices[1] + dimY - 2 == globalPointNum) {
            // std::cout << "CUNTupper" << std::endl;
            #pragma omp parallel for collapse(2)
            for (int i = 1; i < dimX - 1; ++i) {
                for (int k = 1; k < dimZ - 1; ++k) {
                    upperPart1[calcBoundaryIndex(i, k, dimZ)] = getValue(i, dimY - 3, k);
                    // std::cout << upperPart1[calcBoundaryIndex(i, k, dimZ)] << " ";
                }
                // std::cout << std::endl;
            }
        } else {
            #pragma omp parallel for collapse(2)
            for (int i = 1; i < dimX - 1; ++i) {
                for (int k = 1; k < dimZ - 1; ++k) {
                    upperPart1[calcBoundaryIndex(i, k, dimZ)] = getValue(i, dimY - 2, k);
                }
            }
        }

        if (startIndices[0] + dimX - 2 == globalPointNum) {
            // std::cout << "CUNTright" << std::endl;
            #pragma omp parallel for collapse(2)
            for (int j = 1; j < dimY - 1; ++j) {
                for (int k = 1; k < dimZ - 1; ++k) {
                    rightPart1[calcBoundaryIndex(j, k, dimZ)] = getValue(dimX - 3, j, k);
                }
            }
        } else {
            #pragma omp parallel for collapse(2)
            for (int j = 1; j < dimY - 1; ++j) {
                for (int k = 1; k < dimZ - 1; ++k) {
                    rightPart1[calcBoundaryIndex(j, k, dimZ)] = getValue(dimX - 2, j, k);
                }
            }
        }

        if (startIndices[1] == 0) {
            // std::cout << "CUNTlower" << std::endl;
            #pragma omp parallel for collapse(2)
            for (int i = 1; i < dimX - 1; ++i) {
                for (int k = 1; k < dimZ - 1; ++k) {
                    lowerPart1[calcBoundaryIndex(i, k, dimZ)] = getValue(i, 2, k);
                }
            }
        } else {
            #pragma omp parallel for collapse(2)
            for (int i = 1; i < dimX - 1; ++i) {
                for (int k = 1; k < dimZ - 1; ++k) {
                    lowerPart1[calcBoundaryIndex(i, k, dimZ)] = getValue(i, 1, k);
                }
            }
        }

        if (startIndices[0] == 0) {
            // std::cout << "CUNTleft" << std::endl;
            #pragma omp parallel for collapse(2)
            for (int j = 1; j < dimY - 1; ++j) {
                for (int k = 1; k < dimZ - 1; ++k) {
                    leftPart1[calcBoundaryIndex(j, k, dimZ)] = getValue(2, j, k);
                }
            }
        } else {
            #pragma omp parallel for collapse(2)
            for (int j = 1; j < dimY - 1; ++j) {
                for (int k = 1; k < dimZ - 1; ++k) {
                    leftPart1[calcBoundaryIndex(j, k, dimZ)] = getValue(1, j, k);
                }
            }
        }

        #pragma omp parallel for collapse(2)
        for (int i = 1; i < dimX - 1; ++i) {
            for (int j = 1; j < dimY - 1; ++j) {
                frontPart1[calcBoundaryIndex(i, j, dimY)] = getValue(i, j, 1);
            }
        }
        #pragma omp parallel for collapse(2)
        for (int i = 1; i < dimX - 1; ++i) {
            for (int j = 1; j < dimY - 1; ++j) {
                backPart1[calcBoundaryIndex(i, j, dimY)] = getValue(i, j, dimZ - 2);
            }
        }

        // for (int j = 0; j < dimY; ++j) {
        //     for (int k = 0; k < dimZ; ++k) {
        //         std::cout << leftPart1[j * dimZ + k] << " ";
        //     }
        //     std::cout << std::endl;
        // }

        // std::cout << "FUCK" << std::endl;
    }

    void syncBoundary(const std::string &part) {
        if (part == "left") {
            #pragma omp parallel for collapse(2)
            for (int j = 1; j < dimY - 1; ++j) {
                for (int k = 1; k < dimZ - 1; ++k) {
                    setValue(0, j, k, leftPart2[calcBoundaryIndex(j, k, dimZ)]);
                }
            }
        } else if (part == "right") {
            #pragma omp parallel for collapse(2)
            for (int j = 1; j < dimY - 1; ++j) {
                for (int k = 1; k < dimZ - 1; ++k) {
                    setValue(dimX - 1, j, k, rightPart2[calcBoundaryIndex(j, k, dimZ)]);
                }
            }
        } else if (part == "lower") {
            #pragma omp parallel for collapse(2)
            for (int i = 1; i < dimX - 1; ++i) {
                for (int k = 1; k < dimZ - 1; ++k) {
                    setValue(i, 0, k, lowerPart2[calcBoundaryIndex(i, k, dimZ)]);
                }
            }
        } else if (part == "upper") {
            #pragma omp parallel for collapse(2)
            for (int i = 1; i < dimX - 1; ++i) {
                for (int k = 1; k < dimZ - 1; ++k) {
                    setValue(i, dimY - 1, k, upperPart2[calcBoundaryIndex(i, k, dimZ)]);
                }
            }
        } else if (part == "front") {
            #pragma omp parallel for collapse(2)
            for (int i = 1; i < dimX - 1; ++i) {
                for (int j = 1; j < dimY - 1; ++j) {
                    setValue(i, j, 0, frontPart2[calcBoundaryIndex(i, j, dimY)]);
                }
            }
        } else if (part == "back") {
            #pragma omp parallel for collapse(2)
            for (int i = 1; i < dimX - 1; ++i) {
                for (int j = 1; j < dimY - 1; ++j) {
                    setValue(i, j, dimZ - 1, backPart2[calcBoundaryIndex(i, j, dimY)]);
                }
            }
        }
    }

    ~Grid() {
        delete[] values;
        delete[] upperPart1;
        delete[] lowerPart1;
        delete[] rightPart1;
        delete[] leftPart1;
        delete[] backPart1;
        delete[] frontPart1;

        delete[] upperPart2;
        delete[] lowerPart2;
        delete[] rightPart2;
        delete[] leftPart2;
        delete[] backPart2;
        delete[] frontPart2;
    }
};
