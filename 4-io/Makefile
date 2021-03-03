# Compiler
CC = gcc
CFLAGS += -Wall -std=gnu1x
LDLIBS += -lm

ARGS = 1 2 3 4 5
TARGETS = out-cmdline.txt out-stdin.txt out-file.txt

default: $(TARGETS)

out-file.txt: read-file input.txt
	./$< input.txt out-file.txt
	cat $@

out-stdin.txt: read-stdin
	echo $(ARGS) | ./$< > $@
	cat $@

out-cmdline.txt: read-cmdline
	./$< $(ARGS) > $@
	cat $@

read-file: read-file.c input.txt
	$(CC) $< -o $@ $(CFLAGS) $(LDLIBS)
read-stdin: read-stdin.c
read-cmdline: read-cmdline.c

input.txt:
	echo $(ARGS) > $@

.PHONEY:clean
clean:                     # this is "clean" target. it does not depend on anything
	rm *.txt
	find ./ -type f -executable -delete


