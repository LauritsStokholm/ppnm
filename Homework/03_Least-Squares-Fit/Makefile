# Compiler
CC = gcc
CFLAGS += -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS += $$(gsl-config --libs)


#==============================================================================#
# For gnuplot
SCRIPT = log.gpi
FORMAT = svg
IMAGE_VIEWER = ristretto
PLOTTER = gnuplot
void =

#==============================================================================#

# Figure 1
myPlot.$(FORMAT): Rutherford.data myplot.data Makefile
	echo 'set terminal $(FORMAT)' > $(SCRIPT)
	echo 'set output "$@"' >> $(SCRIPT)
	echo 'set key out' >> $(SCRIPT)
	echo 'set title "Ordinary Least-Squares fit by QR-decomposition"' >> $(SCRIPT)
	echo 'set grid' >> $(SCRIPT)
	echo 'set xlabel "Time in days"' >> $(SCRIPT)
	echo 'set ylabel "Activity of ThX (relative units)"' >> $(SCRIPT)
	echo 'plot \$(void)' >> $(SCRIPT)
	echo '"$(word 1, $^)" using 1:2:3 with yerrorbars lc "red" title "Data" \$(void)' >> $(SCRIPT)
	echo ', "$(word 2, $^)" using 1:2 with lines dt 2 lc "black" title "Fmin" \$(void)' >> $(SCRIPT)
	echo ', "$(word 2, $^)" using 1:3 with lines dt 1 lc "black" title "Fval" \$(void)' >> $(SCRIPT)
	echo ', "$(word 2, $^)" using 1:4 with lines dt 2 lc "black" title "Fmax" \$(void)' >> $(SCRIPT)
	$(PLOTTER) $(SCRIPT)
	$(IMAGE_VIEWER) $@
#==============================================================================#

# Output of main
myplot.data: main
	./main
main: main.c


#==============================================================================#
#DATAFILE
#Time t (days)                     : 1,  2,  3, 4, 6, 9,   10,  13,  15
#Activity y of ThX (relative units): 117,100,88,72,53,29.5,25.2,15.2,11.1
Rutherford.data: data
	./$<

data: data.c

#==============================================================================#


.PHONEY: clean
clean:  # this is "clean" target. it does not depend on anything
	find ./ -type f -executable -delete
	rm log.gpi *.data


