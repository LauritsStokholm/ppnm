# Compiler
CC = gcc
CFLAGS += -Wall -pedantic -std=gnu1x $$(gsl-config --cflags)
LDLIBS += $$(gsl-config --libs)

results.txt: main
	./$< > $@
	cat $@

main: main.c

.PHONEY: clean
clean:                     # this is "clean" target. it does not depend on anything
	find ./ -type f -executable -delete
	$(RM) *.txt


