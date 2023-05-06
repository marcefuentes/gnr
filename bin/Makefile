# Compiler and compiler flags
CC = gcc
CFLAGS := -g -std=gnu99 -Wall -Wextra -Wpedantic -Wno-unused-parameter \
          -Wno-unused-function -Wshadow -Wwrite-strings -I../src -I../dtnorm/src

# Optimization and profiling flags
#opt := -DNDEBUG -O3  -finline-functions  # For full optimization
opt := -O0 -fno-inline-functions -DDEBUG     # For debugging
#prof := -pg -rdynamic                    # For profiling
prof :=
CFLAGS += $(opt) $(prof)

# Libraries
LIBS := -L/usr/local/lib -lgsl -lgslcblas -lpthread -lm

# Name of the final executable
TARGET = gnr

# Directories and sources
VPATH = ../src ../dtnorm/src
SOURCES = $(wildcard $(VPATH:%=%/*.c))
OBJS = $(SOURCES:%.c=%.o)
DEPS = $(OBJS:%.o=%.d)

# Generate the GCC includes parameters by adding -I before each source folder
INCLUDES = $(foreach dir, $(VPATH), $(addprefix -I, $(dir)))

# Generate rules to build each object file
%.o: %.c sim.h
	@echo Building $@
	$(CC) $(CFLAGS) $(INCLUDES) -c -o $@ $< -MMD

# Build the final executable
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(LIBS) $(OBJS) -o $(TARGET)

# Include dependencies
-include $(DEPS)

.PHONY: all clean test

all: $(TARGET)

test: $(TARGET)
	-./$(TARGET) t
	@echo Completed all tests.

clean:
	rm -f $(OBJS) $(DEPS) $(TARGET)
	@echo Cleaned.

