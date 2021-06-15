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
ode.$(FORMAT): rk12_ann_ode.txt rk23_ann_ode.txt
	echo 'set terminal eps' > $(SCRIPT)
	echo 'set output "$@"' >> $(SCRIPT)
	echo 'set key out' >> $(SCRIPT)
	echo 'set title "This is a title"' >> $(SCRIPT)
	echo 'set grid' >> $(SCRIPT)
	echo 'set xlabel "x"' >> $(SCRIPT)
	echo 'set ylabel "y"' >> $(SCRIPT)
	echo 'plot \$(void)' >> $(SCRIPT)
	echo '"$(word 1, $^)"   using 1:2  with p title "datapoints rk12 data" \$(void)' >> $(SCRIPT)
	echo ', "$(word 2, $^)" using 1:2  with p title "datapoints rk23 data" \$(void)' >> $(SCRIPT)
	$(PLOTTER) $(SCRIPT)
	gv $@

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


# For comparison in task c
rk12_ann_ode.txt rk23_ann_ode.txt: my_ode
my_ode: my_ode.c

#==============================================================================#

.PHONEY:clean
clean:                     # this is "clean" target. it does not depend on anything
	rm -f $(TARGETS)
	find ./ -type f -executable -delete
	$(RM) *.txt *.eps *.gpi


