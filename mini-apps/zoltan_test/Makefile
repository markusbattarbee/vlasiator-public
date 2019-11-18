
FLAGS = -std=c++11  -g



INCLUDES = -I$(CRAY_TRILINOS_PREFIX_DIR)/include -I$(DCCRG) -L$(CRAY_TRILINOS_PREFIX_DIR)/lib   -lboost_program_options -lzoltan
CXX = CC $(FLAGS)
DCCRG = /proj/vlasov/libraries/dccrg_new_neighbours/

PROGRAMS = load_balancing_test multi_stage_load_balancing

all: $(PROGRAMS)

load_balancing_test: load_balancing_test.cpp $(DCCRG)/dccrg.hpp $(DCCRG)/dccrg_cartesian_geometry.hpp Makefile
	$(CXX) $< $(INCLUDES) -DDEBUG -o $@

d: data
data:
	rm -f *vtk *visit

c: clean
clean: data
	rm -f $(PROGRAMS)

