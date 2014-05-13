program_NAME := shelly
CC = cc
program_C_SRCS := $(wildcard *.c) #$(wildcard util/*.c)
program_C_OBJS := ${program_C_SRCS:.c=.o}
program_OBJS := $(program_C_OBJS)
program_INCLUDE_DIRS := /mt/fielding-home/hemingway/code/fourier-methods/
program_LIBRARY_DIRS := /mt/fielding-home/hemingway/code/fourier-methods/
program_LIBRARIES := pthread gsl gslcblas ffs fftw3 m
#program_FLAGS := -Wall -Wextra -O3 -std=c99 -Wshadow # speed
program_FLAGS := -Wall -Wextra -g -std=c99 -Wshadow # debug

CFLAGS += $(foreach includedir,$(program_INCLUDE_DIRS),-I$(includedir))
CFLAGS += $(program_FLAGS)
LDFLAGS += $(foreach librarydir,$(program_LIBRARY_DIRS),-L$(librarydir))
LDFLAGS += $(foreach library,$(program_LIBRARIES),-l$(library))

.PHONY: all clean distclean

all: $(program_NAME)

# compile with $CC, but link with f95
$(program_NAME): $(program_OBJS)
	cc $(program_OBJS) $(LDFLAGS)  -o $(program_NAME)

clean:
	@- $(RM) $(program_NAME)
	@- $(RM) $(program_OBJS)

distclean: clean
