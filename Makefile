CC = mpicc
CFLAGS = -fopenmp -O3 -Wall
INCLUDES = -Iinclude
LIBS = -lm

SRC = src/stencil_template_parallel.c
TARGET = stencil_hybrid

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $^ $(LIBS)

clean:
	rm -f $(TARGET)
