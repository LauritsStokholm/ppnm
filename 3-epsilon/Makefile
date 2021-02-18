CC := gcc
CFLAGS := -Wall -pedantic -fwrapv -O3

HEADER := main.h
CFILES := main.c
OUTPUT := out.txt main

$(word 1, $(OUTPUT)): $(word 2, $(OUTPUT))
	./$< > $@
	cat $@

$((word 2, $(OUTPUT)): $(HEADER) $(CFILES)

.PHONEY:clean
clean:
	rm $(OUTPUT)

