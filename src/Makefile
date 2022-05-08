CXX := g++-11
CXXFLAGS = -std=gnu++17 -Wall -Wextra -O3

SRCS   = $(wildcard *.cpp)
OBJS   = $(SRCS:.cpp=.o)

main: $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

%: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

graph.o: graph.hpp

clean:
	rm *.o main

.PHONY: clean