#include <string>
#include <chrono>
#include <iostream>
#include <cmath>
#include <mpi.h>

#include "grid.h"
#include "grid_updater.h"

#define MAIN_PROC 0

void sendBoundaries(MPI_Comm gridComm, Grid &prevGrid) {
    MPI_Status status[12];
    MPI_Request requests[12];
    int neighbors[6];

    MPI_Cart_shift(gridComm, 0, 1, &neighbors[0], &neighbors[1]); // left, right
    MPI_Cart_shift(gridComm, 1, 1, &neighbors[2], &neighbors[3]); // lower, upper
    MPI_Cart_shift(gridComm, 2, 1, &neighbors[4], &neighbors[5]); // front, back

    int requestIndex = 0;

    if (neighbors[0] != MPI_PROC_NULL) {
        MPI_Isend(prevGrid.leftPart1, prevGrid.dimY * prevGrid.dimZ, MPI_DOUBLE, neighbors[0], 0, gridComm, &requests[requestIndex++]);
        MPI_Irecv(prevGrid.leftPart2, prevGrid.dimY * prevGrid.dimZ, MPI_DOUBLE, neighbors[0], 1, gridComm, &requests[requestIndex++]);
    }

    if (neighbors[1] != MPI_PROC_NULL) {
        MPI_Isend(prevGrid.rightPart1, prevGrid.dimY * prevGrid.dimZ, MPI_DOUBLE, neighbors[1], 1, gridComm, &requests[requestIndex++]);
        MPI_Irecv(prevGrid.rightPart2, prevGrid.dimY * prevGrid.dimZ, MPI_DOUBLE, neighbors[1], 0, gridComm, &requests[requestIndex++]);
    }

    if (neighbors[2] != MPI_PROC_NULL) {
        MPI_Isend(prevGrid.lowerPart1, prevGrid.dimX * prevGrid.dimZ, MPI_DOUBLE, neighbors[2], 3, gridComm, &requests[requestIndex++]);
        MPI_Irecv(prevGrid.lowerPart2, prevGrid.dimX * prevGrid.dimZ, MPI_DOUBLE, neighbors[2], 2, gridComm, &requests[requestIndex++]);
    }

    if (neighbors[3] != MPI_PROC_NULL) {
        MPI_Isend(prevGrid.upperPart1, prevGrid.dimX * prevGrid.dimZ, MPI_DOUBLE, neighbors[3], 2, gridComm, &requests[requestIndex++]);
        MPI_Irecv(prevGrid.upperPart2, prevGrid.dimX * prevGrid.dimZ, MPI_DOUBLE, neighbors[3], 3, gridComm, &requests[requestIndex++]);
    }

    if (neighbors[4] != MPI_PROC_NULL) {
        MPI_Isend(prevGrid.frontPart1, prevGrid.dimX * prevGrid.dimY, MPI_DOUBLE, neighbors[4], 5, gridComm, &requests[requestIndex++]);
        MPI_Irecv(prevGrid.frontPart2, prevGrid.dimX * prevGrid.dimY, MPI_DOUBLE, neighbors[4], 4, gridComm, &requests[requestIndex++]);
    }

    if (neighbors[5] != MPI_PROC_NULL) {
        MPI_Isend(prevGrid.backPart1, prevGrid.dimX * prevGrid.dimY, MPI_DOUBLE, neighbors[5], 4, gridComm, &requests[requestIndex++]);
        MPI_Irecv(prevGrid.backPart2, prevGrid.dimX * prevGrid.dimY, MPI_DOUBLE, neighbors[5], 5, gridComm, &requests[requestIndex++]);
    }

    MPI_Waitall(requestIndex, requests, status);

    if (neighbors[0] != MPI_PROC_NULL) {
        prevGrid.syncBoundary("left");
    }

    if (neighbors[1] != MPI_PROC_NULL) {
        prevGrid.syncBoundary("right");
    }

    if (neighbors[2] != MPI_PROC_NULL) {
        prevGrid.syncBoundary("lower");
    }

    if (neighbors[3] != MPI_PROC_NULL) {
        prevGrid.syncBoundary("upper");
    }

    if (neighbors[4] != MPI_PROC_NULL) {
        prevGrid.syncBoundary("front");
    }

    if (neighbors[5] != MPI_PROC_NULL) {
        prevGrid.syncBoundary("back");
    }
}

double parseLengthArgument(const std::string& arg) {
    if (arg == "pi") {
        return M_PI;
    }
    return std::stod(arg);
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int procNum, procID;
    MPI_Comm_size(MPI_COMM_WORLD, &procNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &procID);

    size_t spatialPointNum = std::stoi(argv[1]);
    size_t timePointNum = std::stoi(argv[2]);
    double Lx = parseLengthArgument(argv[3]);
    double Ly = parseLengthArgument(argv[4]);
    double Lz = parseLengthArgument(argv[5]);
    double T = std::stod(argv[6]);

    int procDims[3] = { 0, 0, 0 };
    MPI_Dims_create(procNum, 3, procDims);

    MPI_Comm gridComm;
    int periods[3] = { 1, 1, 0 };
    MPI_Cart_create(MPI_COMM_WORLD, 3, procDims, periods, false, &gridComm);

    int procCoords[3];
    MPI_Cart_coords(gridComm, procID, 3, procCoords);

    size_t blockStartIndices[3], blockEndIndices[3], blockSize[3];
    for (int dim = 0; dim < 3; ++dim) {
        size_t baseBlockSize = spatialPointNum / procDims[dim];
        size_t remainder = spatialPointNum % procDims[dim];
        
        blockStartIndices[dim] = baseBlockSize * procCoords[dim];
        blockEndIndices[dim] = baseBlockSize * (procCoords[dim] + 1);

        if (procCoords[dim] == procDims[dim] - 1) {
            blockEndIndices[dim] += remainder;
        }

        blockSize[dim] = blockEndIndices[dim] - blockStartIndices[dim];
        blockSize[dim] += 2 ;
    }

    double startTime = MPI_Wtime();

    double timeStep = T / (double) timePointNum;

    GridUpdater gridUpdater(Lx, Ly, Lz, timeStep, spatialPointNum);

    Grid prevGrid = Grid(spatialPointNum, Lx, Ly, Lz, blockSize, blockStartIndices);
    Grid curGrid = Grid(spatialPointNum, Lx, Ly, Lz, blockSize, blockStartIndices);
    Grid nextGrid = Grid(spatialPointNum, Lx, Ly, Lz, blockSize, blockStartIndices);

    gridUpdater.initializeGridT0(prevGrid, blockStartIndices);
    prevGrid.fillBoundary();
    sendBoundaries(gridComm, prevGrid);

    MPI_Barrier(MPI_COMM_WORLD);

    double maxError = gridUpdater.initializeGridT1(prevGrid, curGrid, blockStartIndices);

    double totalMaxError;
    MPI_Reduce(&maxError, &totalMaxError, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (procID == MAIN_PROC) {
        std::cout << "Макс. ошибка: " << totalMaxError << std::endl;
    }

    double curTotalMaxError = 0.0;
    for (size_t t = 2; t < timePointNum; ++t) {
        curGrid.fillBoundary();
        sendBoundaries(gridComm, curGrid);

        double maxError = gridUpdater.updateGrid(prevGrid, curGrid, nextGrid, t * timeStep, blockStartIndices);
        curTotalMaxError = std::max(curTotalMaxError, maxError);

        MPI_Reduce(&maxError, &totalMaxError, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if (procID == MAIN_PROC) {
            std::cout << t << " шаг:  Макс. ошибка: " << totalMaxError << std::endl;
        }
        
        gridUpdater.move(prevGrid, curGrid);
        gridUpdater.move(curGrid, nextGrid);
    }

    double endTime = MPI_Wtime();
    double procWorkTime = endTime - startTime;
    double maxProcWorkTime;
    MPI_Reduce(&procWorkTime, &maxProcWorkTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    MPI_Reduce(&curTotalMaxError, &totalMaxError, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (procID == MAIN_PROC) {
        std::cout << "Расчет завершен.\n";
        std::cout << "Время выполнения: " << maxProcWorkTime * 1000  << " ms.\n";
        std::cout << "Максимальная ошибка на сетке: " << totalMaxError << "\n";
    }

    MPI_Finalize();

    return 0;
}