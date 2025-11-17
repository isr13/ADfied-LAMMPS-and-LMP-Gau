DCO_BASE_DIR=../..
CXXFLAGS=-Wall -pedantic -Wextra -DNDEBUG -Werror -O3 -I$(DCO_BASE_DIR)/include -DDCO_DISABLE_AUTO_WARNING -DDCO_DISABLE_AVX2_WARNING -fopenmp -march=native 
CXX=g++

source=main.cpp
objects=$(source:.cpp=.o)

all : main.exe run check 

main.exe : $(objects)
	$(CXX) $^  $(DCO_BASE_DIR)/lib/libdcoc.a -ldl -o main.exe -fopenmp

run : main.exe
	./main.exe  > out

check : run
	-diff --strip-trailing-cr out ref_results

generate_results : main.exe
	./main.exe  > ref_results

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $^

clean:
	rm -f *.o main.exe out
