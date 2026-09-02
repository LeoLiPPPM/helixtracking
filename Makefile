CXX := g++
CXXFLAGS := -std=c++20 -O2 -Wall -Wextra -Wpedantic -Wconversion -Iinclude
BUILD_DIR := build
BIN_DIR := bin

.PHONY: all test demo clean

all: $(BIN_DIR)/helixtracking

$(BUILD_DIR)/track.o: src/track.cpp include/helixtracking/track.hpp
	@mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BIN_DIR)/helixtracking: src/main.cpp $(BUILD_DIR)/track.o
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(BIN_DIR)/test_track: tests/test_track.cpp $(BUILD_DIR)/track.o
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $^ -o $@

test: $(BIN_DIR)/test_track
	./$(BIN_DIR)/test_track

demo: $(BIN_DIR)/helixtracking
	./$(BIN_DIR)/helixtracking --output outputs --seed 20260902

clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)

