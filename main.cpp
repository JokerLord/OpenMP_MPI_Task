#include <string>
#include <chrono>
#include <iostream>
#include <cmath>

#include "grid.h"
#include "grid_updater.h"


int main(int argc, char *argv[]) {
    size_t spatialPointNum = std::stoi(argv[1]);
    size_t timePointNum = std::stoi(argv[2]);
    double Lx = std::stod(argv[3]);
    double Ly = std::stod(argv[4]);
    double Lz = std::stod(argv[5]);
    double T = std::stod(argv[6]);

    std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();

    double timeStep = T / (double) timePointNum;

    GridUpdater gridUpdater(Lx, Ly, Lz, timeStep);

    Grid prevGrid = Grid(spatialPointNum, Lx, Ly, Lz);
    Grid curGrid = Grid(spatialPointNum, Lx, Ly, Lz);
    Grid nextGrid = Grid(spatialPointNum, Lx, Ly, Lz);

    gridUpdater.initializeGridT0(prevGrid);
    double maxError = gridUpdater.initializeGridT1(prevGrid, curGrid);

    std::cout << "Макс. ошибка: " << maxError << std::endl;

    double maxErrorTotal = 0.0;

    for (size_t t = 2; t < timePointNum; ++t) {
        double maxError = gridUpdater.updateGrid(prevGrid, curGrid, nextGrid, t * timeStep);
        maxErrorTotal = std::max(maxErrorTotal, maxError);
        std::cout << "Макс. ошибка: " << maxError << std::endl;

        gridUpdater.move(prevGrid, curGrid);
        gridUpdater.move(curGrid, nextGrid);
    }

    std::chrono::steady_clock::time_point endTime = std::chrono::steady_clock::now();
    std::chrono::duration<double, std::milli> workTime = endTime - startTime;

    std::cout << "Расчет завершен.\n";
    std::cout << "Время выполнения: " << workTime.count() << "ms.\n";
    std::cout << "Максимальная ошибка на сетке: " << maxErrorTotal << "\n";

    return 0;
}