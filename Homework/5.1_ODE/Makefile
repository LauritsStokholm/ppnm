# Compiler
CC = gcc
CFLAGS += -Wall -pedantic -std=gnu1x $$(gsl-config --cflags)
LDLIBS += $$(gsl-config --libs)
#==============================================================================#
# For gnuplot
SCRIPT = log.gpi
FORMAT = eps
PLOTTER = gnuplot
void =

#==============================================================================#
harmonic.$(FORMAT): task_a.txt Makefile
	echo 'set terminal eps' > $(SCRIPT)
	echo 'set output "$@"' >> $(SCRIPT)
	echo 'set key out' >> $(SCRIPT)
	echo 'set title "This is a title"' >> $(SCRIPT)
	echo 'set grid' >> $(SCRIPT)
	echo 'set xlabel "x"' >> $(SCRIPT)
	echo 'set ylabel "y"' >> $(SCRIPT)
	echo 'plot \$(void)' >> $(SCRIPT)
	echo '"$(word 1, $^)"   using 1:2 with l dt 1 lc "black" title "rk12, y0" \$(void)' >> $(SCRIPT)
	echo ', "$(word 1, $^)" using 1:3 with l dt 1 lc "black" title "rk12, y1" \$(void)' >> $(SCRIPT)
	echo ', "$(word 1, $^)" using 1:4 with l dt 2 lc "red" title "rk12, cos" \$(void)' >> $(SCRIPT)
	echo ', "$(word 1, $^)" using 1:5 with l dt 2 lc "blue" title "rk12, sin" \$(void)' >> $(SCRIPT)
	$(PLOTTER) $(SCRIPT)
	gv $@

task_a.txt: main
	./$<

main: main.c

.PHONEY: clean
clean:
	rm main *.txt *.eps

