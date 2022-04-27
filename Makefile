CC     := g++-11
CFLAGS := -O3 -std=c++20 -I./include -fopenmp -Wall
SRCDIR := src
OBJDIR := obj
SRCS   := $(wildcard $(SRCDIR)/*.cpp)
OBJS   := $(subst $(SRCDIR)/, $(OBJDIR)/, $(patsubst %.cpp,%.o,$(SRCS)))
TARGET := a.out

$(TARGET): $(OBJS) main.cpp
	$(CC) $(CFLAGS) $^ -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CC) $(CFLAGS) -c $^ -o $@
	
.PHONY: clean
clean:
	rm -f ./obj/*
