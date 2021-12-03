CXX 			:= mpic++
CXX_FLAGS 		:= -std=c++2a -Wall
BIN				:= bin
SRC				:= src

INC_DIRS		:= /usr/include/hdf5/serial include
INCLUDE 		:= $(foreach d, $(INC_DIRS), -I$d)

LIB_DIRS		:= /usr/lib/x86_64-linux-gnu/hdf5/serial/lib lib
LIB		 		:= $(foreach d, $(LIB_DIRS), -L$d)

LIB_FLAGS		:= -lhdf5

EXECUTABLE		:= mld

all: $(BIN)/$(EXECUTABLE)

$(BIN)/$(EXECUTABLE): $(SRC)/*.cpp
	@echo "Building..."
	@mkdir -p $(BIN)
	$(CXX) $(CXX_FLAGS) $(INCLUDE) $(LIB) $^ -o $@ $(LIB_FLAGS)
