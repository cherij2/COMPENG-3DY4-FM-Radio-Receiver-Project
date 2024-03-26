# this makefile is intended for g++ on Linux

CC = g++
CFLAGS = -c -Wall -O3 -std=c++17 #-mcpu=cortex-a72+crypto #comment -mcpu=cortex-a72+crypto if running on VM
LDFLAGS = -pthread
INCLUDE = -I../include/
SOURCES = project.cpp iofunc.cpp filter.cpp  fourier.cpp  genfunc.cpp  logfunc.cpp RFfront.cpp mode.cpp thread.cpp
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = project

all: $(EXECUTABLE) clean

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDE) $< -o $@

clean:
	-rm $(OBJECTS)
