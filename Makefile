# Имя исполняемого файла
TARGET = main.out

# Компилятор и флаги
CXX = mpixlC
CXXFLAGS = -std=c++11 -O3 -fopenmp

# Список файлов
SRC = main.cpp
HEADERS = grid.h grid_updater.h

# Правило по умолчанию для сборки исполняемого файла
all: $(TARGET)

# Правило для сборки исполняемого файла
$(TARGET): $(SRC) $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)

# Очистка сборки
clean:
	rm -f $(TARGET)
