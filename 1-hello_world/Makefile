# Compiler
CC = gcc
CFLAGS += -Wall -std=gnu1x# $$(gsl-config --cflags)
LDLIBS += -lm #$$(gsl-config --libs)

# TARGET, EXECUTABLE
TAR = out.txt
EXE = hello

#==============================================================================
out.txt: hello
	./$< > $@
	cat $@

hello: hello.c

.PHONEY: clean
clean:
	rm $(TAR) $(EXE)
