CC = g++
INCDIR = ../include
CFLAGS = -O3 -Wall -std=c++11 

TPSA_LIBS = tpsaso
TRACKING_LIBS = map_method
LIBS = -L. -l$(TRACKING_LIBS) -lstdc++ -lm

SRC = $(wildcard *.cc)
TARGET = $(SRC:.cc=)
DEPS = $(wildcard $(INCDIR)/*.hpp)
DEPS += $(wildcard $(INCDIR)/*.h)

.PHONY: all
all: $(TARGET)

.cc:
	$(CC) -o $@ $< $(CFLAGS) $(LIBS)

.PHONY: clean
clean:
	rm -f $(TARGET)