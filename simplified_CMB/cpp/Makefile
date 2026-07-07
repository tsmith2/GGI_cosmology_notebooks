CXX ?= c++
CXXFLAGS ?= -O3 -DNDEBUG -std=c++17 -Wall -Wextra -pedantic

TARGET = two_fluid_tt
SRC = two_fluid_tt.cpp

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $@ $<

clean:
	rm -f $(TARGET)
