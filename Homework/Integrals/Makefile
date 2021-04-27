# Compiler
CC = gcc
CFLAGS += -Wall -pedantic -std=gnu1x $$(gsl-config --cflags)
LDLIBS += $$(gsl-config --libs)

# For big projects:
# VPATH is included for all files
#VPATH = ~~@~~
## vpath is file specific
#vpath .mca ..~//~


#==============================================================================#
# For gnuplot
SCRIPT = log.gpi
FORMAT = eps
PLOTTER = gnuplot
void =

#==============================================================================#
TARGETS = myPlot.eps output.txt
# Remember: $^ (prerequisite) and $@ (target)
#==============================================================================#

#.PHONEY: all
all: $(TARGETS)
#==============================================================================#

# Figure 1
myPlot.$(FORMAT): output.txt Makefile
	echo 'set terminal eps' > $(SCRIPT)
	echo 'set output "$@"' >> $(SCRIPT)
	echo 'set key out' >> $(SCRIPT)
	echo 'set title "This is a title"' >> $(SCRIPT)
	echo 'set grid' >> $(SCRIPT)
	echo 'set xlabel "x"' >> $(SCRIPT)
	echo 'set ylabel "y"' >> $(SCRIPT)
	echo 'plot \$(void)' >> $(SCRIPT)
	echo '"$(word 1, $^)" using 1:2 with points title "Data" \$(void)' >> $(SCRIPT)
	echo ', "$(word 1, $^)" using 1:3 with lines title "fit"' >> $(SCRIPT)
	$(PLOTTER) $(SCRIPT)
	gv $@
#==============================================================================#

# Output of main
main: main.c
#==============================================================================#

.PHONEY:clean
clean:                     # this is "clean" target. it does not depend on anything
	rm -f $(TARGETS)
	find ./ -type f -executable -delete
	$(RM) *.dat $(SCRIPT) *.svg *.png *.pdf log* *.eps


