INCLUDES=
OBJECTS=environment.o
TARGET=environment_packed
HYPRE_LIB=-L/opt/hypre/2.9.0b/lib -lHYPRE -lblas -llapack
HYPRE_INC=-I/opt/hypre/2.9.0b/include
default: $(TARGET)

%.o: %.c $(HEADERS)
	mpicc $(HYPRE_INC) -c $< -o $@

$(TARGET): $(OBJECTS)
	mpicc -fopenmp $(OBJECTS) $(HYPRE_LIB) -lm -o $@

clean:
	rm -f $(OBJECTS)
	rm -f $(TARGET)
