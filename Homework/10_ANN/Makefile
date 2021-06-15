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
# Figure task a
fit.$(FORMAT): datapoints.txt fit.txt Makefile
	echo 'set terminal eps' > $(SCRIPT)
	echo 'set output "$@"' >> $(SCRIPT)
	echo 'set key out' >> $(SCRIPT)
	echo 'set title "This is a title"' >> $(SCRIPT)
	echo 'set grid' >> $(SCRIPT)
	echo 'set xlabel "x"' >> $(SCRIPT)
	echo 'set ylabel "y"' >> $(SCRIPT)
	echo 'plot \$(void)' >> $(SCRIPT)
	echo '"$(word 1, $^)"   using 1:2  with p title "datapoints f" \$(void)' >> $(SCRIPT)
	echo ', "$(word 1, $^)" using 1:3  with p title "datapoints df" \$(void)' >> $(SCRIPT)
	echo ', "$(word 1, $^)" using 1:4  with p title "datapoints F" \$(void)' >> $(SCRIPT)
	echo ', "$(word 2, $^)" using 1:2  with l title "ann fitted f" \$(void)' >> $(SCRIPT)
	echo ', "$(word 2, $^)" using 1:3  with l title "ann fitted df" \$(void)' >> $(SCRIPT)
	echo ', "$(word 2, $^)" using 1:4  with l title "ann fitted F" \$(void)' >> $(SCRIPT)
	$(PLOTTER) $(SCRIPT)
	gv $@

#==============================================================================#
# Output of main
datapoints.txt fit.txt: main
	./main

main: main.c
#==============================================================================#

.PHONEY:clean
clean:                     # this is "clean" target. it does not depend on anything
	rm -f $(TARGETS)
	find ./ -type f -executable -delete
	$(RM) *.txt *.eps *.gpi


