# Compiler
CC = gcc
CFLAGS += -Wall -std=gnu1x # $$(gsl-config --cflags)
LDLIBS += -lm #$$(gsl-config --libs)

out.txt: main
	./$< > $@
	cat $@

main: main.o komplex.o
main.o komplex.o: komplex.h

.PHONEY: clean
clean:
	$(RM) main main.o komplex.o out.txt
