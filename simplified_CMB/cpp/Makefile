CXX ?= c++
CXXFLAGS ?= -O3 -DNDEBUG -std=c++17 -Wall -Wextra -pedantic

TARGET = two_fluid_tt
SRC = two_fluid_tt.cpp
BESSEL_CACHE_PATTERNS = *.bessel *.besselx *bessel*.bin *bessel*.cache

.PHONY: all clean rebuild

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $@ $<

clean:
	rm -f $(TARGET)
	rm -f $(BESSEL_CACHE_PATTERNS)
	rm -f $(addprefix ../.cache/,$(BESSEL_CACHE_PATTERNS))

rebuild: clean all
