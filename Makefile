CXX := g++
CXXFLAGS := -Wall -Wextra -O2 -std=c++17

TARGET := shade

SRCS := proj2.cpp
OBJS := $(SRCS:.cpp=.o)

all: $(TARGET)


$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)


%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@


clean:
	rm -f $(OBJS) $(TARGET)

