# Compiler
CC = gcc
CFLAGS += -Wall -pedantic -std=gnu1x $$(gsl-config --cflags)
LDLIBS += $$(gsl-config --libs)

#==============================================================================#
# For gnuplot
SCRIPT = log.gpi
FORMAT = svg
IMAGE_VIEWER = ristretto
PLOTTER = gnuplot
void =
#==============================================================================#
# Default
all: threebody.$(FORMAT) SIR.$(FORMAT) harmonic.$(FORMAT)


# Figure task a
threebody.$(FORMAT): rk12_threebody.txt rk23_threebody.txt Makefile
	echo 'set terminal $(FORMAT)' > $(SCRIPT)
	echo 'set output "$@"' >> $(SCRIPT)
	echo 'set key out' >> $(SCRIPT)
	echo 'set title "This is a title"' >> $(SCRIPT)
	echo 'set grid' >> $(SCRIPT)
	echo 'set xlabel "x"' >> $(SCRIPT)
	echo 'set ylabel "y"' >> $(SCRIPT)
	echo 'plot \$(void)' >> $(SCRIPT)
	echo '"$(word 1, $^)"   using 2:3  with lp title "rk12, planet 1" \$(void)' >> $(SCRIPT)
	echo ', "$(word 1, $^)" using 4:5  with lp title "rk12, planet 2" \$(void)' >> $(SCRIPT)
	echo ', "$(word 1, $^)" using 6:7  with lp title "rk12, planet 3" \$(void)' >> $(SCRIPT)
	echo ', "$(word 2, $^)" using 2:3  with lp title "rk23, planet 1" \$(void)' >> $(SCRIPT)
	echo ', "$(word 2, $^)" using 4:5  with lp title "rk23, planet 2" \$(void)' >> $(SCRIPT)
	echo ', "$(word 2, $^)" using 6:7  with lp title "rk23, planet 3" \$(void)' >> $(SCRIPT)
	$(PLOTTER) $(SCRIPT)
	$(IMAGE_VIEWER) $@

#==============================================================================#
# Figure task a
SIR.$(FORMAT): rk12_SIR.txt rk23_SIR.txt Makefile
	echo 'set terminal $(FORMAT)' > $(SCRIPT)
	echo 'set output "$@"' >> $(SCRIPT)
	echo 'set key out' >> $(SCRIPT)
	echo 'set title "This is a title"' >> $(SCRIPT)
	echo 'set grid' >> $(SCRIPT)
	echo 'set xlabel "x"' >> $(SCRIPT)
	echo 'set ylabel "y"' >> $(SCRIPT)
	echo 'plot \$(void)' >> $(SCRIPT)
	echo '"$(word 1, $^)" using 1:2 with lp title "rk12, S" \$(void)' >> $(SCRIPT)
	echo ', "$(word 1, $^)" using 1:3 with lp title "rk12, I" \$(void)' >> $(SCRIPT)
	echo ', "$(word 1, $^)" using 1:4 with lp title "rk12, R" \$(void)' >> $(SCRIPT)
	echo ', "$(word 2, $^)" using 1:2 with lp title "rk23, S" \$(void)' >> $(SCRIPT)
	echo ', "$(word 2, $^)" using 1:3 with lp title "rk23, I" \$(void)' >> $(SCRIPT)
	echo ', "$(word 2, $^)" using 1:4 with lp title "rk23, R" \$(void)' >> $(SCRIPT)
	$(PLOTTER) $(SCRIPT)
	$(IMAGE_VIEWER) $@


harmonic.$(FORMAT): rk12_harmonic.txt rk23_harmonic.txt Makefile
	echo 'set terminal $(FORMAT)' > $(SCRIPT)
	echo 'set output "$@"' >> $(SCRIPT)
	echo 'set key out' >> $(SCRIPT)
	echo 'set title "This is a title"' >> $(SCRIPT)
	echo 'set grid' >> $(SCRIPT)
	echo 'set xlabel "x"' >> $(SCRIPT)
	echo 'set ylabel "y"' >> $(SCRIPT)
	echo 'plot \$(void)' >> $(SCRIPT)
	echo '"$(word 1, $^)" using 1:2 with lp title "rk12, y0" \$(void)' >> $(SCRIPT)
	echo ', "$(word 1, $^)" using 1:3 with lp title "rk12, y1" \$(void)' >> $(SCRIPT)
	echo ', "$(word 2, $^)" using 1:2 with lp title "rk23, y0" \$(void)' >> $(SCRIPT)
	echo ', "$(word 2, $^)" using 1:3 with lp title "rk23, y1" \$(void)' >> $(SCRIPT)
	$(PLOTTER) $(SCRIPT)
	$(IMAGE_VIEWER) $@

#==============================================================================#

# Output of main

task_a.txt rk12_harmonic.txt rk23_harmonic.txt rk12_SIR.txt rk23_SIR.txt rk12_threebody.txt rk23_threebody.txt: main
	./main

main: main.c
#==============================================================================#

.PHONEY:clean
clean:                     # this is "clean" target. it does not depend on anything
	rm -f $(TARGETS)
	find ./ -type f -executable -delete
	$(RM) *.dat *.txt $(SCRIPT) *.svg *.png *.pdf log* *.eps


