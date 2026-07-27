# Generalized Makefile for C++ project
# Compiler: g++-13
# Note: Upgraded to C++17 because modern googletest requires it

CXX      = g++-13
CXXFLAGS = -std=c++17 -Wall -Wextra -O2 -g
INCLUDES = -I./include -I./libraries/eigen -I./libraries/mavlink
LDFLAGS  =

# ====================== Main Application ======================
SRCDIR   = source
OBJDIR   = build
SOURCES  = $(wildcard $(SRCDIR)/*.cpp)
OBJECTS  = $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SOURCES))
TARGET   = spacecraft_simulator

# Identify main object so we can exclude it from tests
MAIN_OBJ := $(OBJDIR)/main.o
APP_OBJECTS := $(filter-out $(MAIN_OBJ), $(OBJECTS))

all: $(TARGET)

$(OBJDIR):
	mkdir -p $(OBJDIR)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp | $(OBJDIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

$(TARGET): $(OBJECTS)
	$(CXX) $(OBJECTS) $(LDFLAGS) -o $@

# ====================== Google Test ======================
GTEST_DIR     := libraries/googletest/googletest
GTEST_INCLUDE := -I$(GTEST_DIR)/include
GTEST_SRC     := $(GTEST_DIR)/src

# Pre-compile googletest (correct flags — no extra -I src)
$(OBJDIR)/gtest-all.o: $(GTEST_SRC)/gtest-all.cc | $(OBJDIR)
	$(CXX) $(CXXFLAGS) $(GTEST_INCLUDE) -I$(GTEST_DIR) -c $< -o $@

$(OBJDIR)/gtest_main.o: $(GTEST_SRC)/gtest_main.cc | $(OBJDIR)
	$(CXX) $(CXXFLAGS) $(GTEST_INCLUDE) -I$(GTEST_DIR) -c $< -o $@

# Test sources (put your test files in tests/)
TEST_SRCDIR  := tests
TEST_SOURCES := $(wildcard $(TEST_SRCDIR)/*.cpp)
TEST_OBJECTS := $(patsubst $(TEST_SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(TEST_SOURCES))

# Compile test files
$(OBJDIR)/%.o: $(TEST_SRCDIR)/%.cpp | $(OBJDIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(GTEST_INCLUDE) -c $< -o $@

TEST_TARGET := test_spacecraft_simulator
GTEST_FLAGS := --gtest_color=yes --gtest_print_time=0

# Build test binary — automatically includes ALL your source objects except main.o
$(TEST_TARGET): $(TEST_OBJECTS) $(OBJDIR)/gtest-all.o $(OBJDIR)/gtest_main.o $(APP_OBJECTS)
	$(CXX) $(CXXFLAGS) $(GTEST_INCLUDE) $^ $(LDFLAGS) -o $@ -pthread

# Run tests
test: $(TEST_TARGET)
	./$(TEST_TARGET) $(GTEST_FLAGS)

# ====================== Housekeeping ======================
clean:
	rm -rf $(OBJDIR) $(TARGET) $(TEST_TARGET) sim.csv debug.txt plots/

run: $(TARGET)
	./$(TARGET)

.PHONY: all clean run test
