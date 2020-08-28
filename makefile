EXEC 		= extract_shapelets
CC			= g++
CFLAGS 		= -lm -Wall -O3 
SRC_DIR 	= ./src
BUILD_DIR 	= ./build
BIN_DIR 	= ./bin
SOURCES		= shapelet_transform.cpp extract_shapelets.cpp

# create the obj variable by substituting the extension of the sources
# and adding a path
_OBJ = $(SOURCES:.cpp=.o) #changes the .cpp extension to .o extension from source file list
OBJ = $(patsubst %,$(BUILD_DIR)/%,$(_OBJ))

all: $(BIN_DIR)/$(EXEC)

$(BIN_DIR)/$(EXEC): $(OBJ)
	$(CC)  -o $@ $^ $(CFLAGS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp 
	$(CC)  -c -o $@ $< $(CFLAGS)

.PHONY: clean
clean:
	$(RM) *.o $(BIN_DIR)/$(EXEC) $(OBJ)
