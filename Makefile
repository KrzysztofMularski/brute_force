CC = clang++ -std=c++17
CFLAGS = -g -Wall -Wextra -Ofast
# LIBS = -DEIGEN_DONT_PARALLELIZE -fopenmp=libomp
LIBS = -fopenmp=libomp

TARGET = brute_force

.PHONY: $(TARGET)
all: $(TARGET)

$(TARGET): $(TARGET).cpp
	$(CC) $(CFLAGS) $(LIBS) $(TARGET).cpp -o $(TARGET)

clean:
	$(RM) $(TARGET)

