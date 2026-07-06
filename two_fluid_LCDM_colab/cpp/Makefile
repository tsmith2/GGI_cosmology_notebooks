CXX ?= c++
CXXFLAGS ?= -O3 -DNDEBUG -std=c++17 -Wall -Wextra -pedantic

# Optional local-machine optimizations:
#   make NATIVE=1        adds -march=native
#   make LTO=1           adds -flto
#
# -pthread is intentionally not enabled by default: this solver is
# single-threaded, while the MCMC parallelism is handled by Python processes.
ifeq ($(NATIVE),1)
CXXFLAGS += -march=native
endif

ifeq ($(LTO),1)
CXXFLAGS += -flto
endif

TARGET = two_fluid_tt
SRC = two_fluid_tt.cpp

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $@ $<

clean:
	rm -f $(TARGET)
