CC := gcc
CFLAGS += -Wall -pedantic -std=gnu1x $$(gsl-config --cflags)
LDLIBS += $$(gsl-config --libs)
EXEC = main
TARGETS = output_linalg.txt output_eigenv.txt

default: output_linalg.txt output_eigenv.txt
	cat $<
	cat $(word 2, $^)


output_eigenv.txt: main
	./$<

output_linalg.txt: main input_linalg.txt
	./$<

main: main.c

input_linalg.txt:
	echo "6.13 -2.90 5.86 8.08 -6.31 -3.89 -4.36 1.00 0.19" > $@


.PHONEY: clean
clean:
	rm *.txt $(EXEC)

