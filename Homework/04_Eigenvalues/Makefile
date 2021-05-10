# Compiler
CC = gcc
CFLAGS += -Wall -pedantic -fwrapv -O3 -std=gnu1x $$(gsl-config --cflags)
LDLIBS += $$(gsl-config --libs)

#==============================================================================#
# For gnuplot
SCRIPT = log.gpi
FORMAT = eps
PLOTTER = gnuplot
void =
#==============================================================================#

# CONVERGENCE TESTS:
# Figure
operations.$(FORMAT): convergence_input.txt convergence_gsl_input.txt Makefile
	echo 'set terminal eps' > $(SCRIPT)
	echo 'set output "$@"' >> $(SCRIPT)
	echo 'set key left' >> $(SCRIPT)
	echo 'set title "This is a title"' >> $(SCRIPT)
	echo 'set grid' >> $(SCRIPT)
	echo 'set xlabel "x"' >> $(SCRIPT)
	echo 'set ylabel "y"' >> $(SCRIPT)
	echo 'plot \$(void)' >> $(SCRIPT)
	echo '"$(word 1, $^)"   u 1:2 w p pt 1 ps 1 t "My implementation"  \$(void)' >> $(SCRIPT)
	echo ', "$(word 2, $^)" u 1:2 w p pt 2 ps 1 t "GSL"\$(void)' >> $(SCRIPT)
	$(PLOTTER) $(SCRIPT)
	gv $@



convergence_gsl_input.txt: gsl_timer
	for i in $$(seq 20 5 200);\
		do\
		time --format "$$i\t%R" --append --output $@ ./$< -n $$i;\
		done

convergence_input.txt: my_timer
	for i in $$(seq 20 5 200);\
		do\
		time --format "$$i\t%R" --append --output $@ ./$< -n $$i;\
		done


# Stripped down for timing only part a
gsl_timer: gsl_timer.c
my_timer: my_timer.c


#==============================================================================#

# Figure 1
myPlot.$(FORMAT): plot.txt Makefile
	echo 'set terminal eps' > $(SCRIPT)
	echo 'set output "$@"' >> $(SCRIPT)
	echo 'set key out' >> $(SCRIPT)
	echo 'set title "Standing Waves with BDC"' >> $(SCRIPT)
	echo 'set grid' >> $(SCRIPT)
	echo 'set xlabel "x"' >> $(SCRIPT)
	echo 'set ylabel "y"' >> $(SCRIPT)
	echo 'plot \$(void)' >> $(SCRIPT)
	echo '"$(word 1, $^)"   u 1:2 w lp ps 0.5 lc "black" dt 1 t "1st" \$(void)' >> $(SCRIPT)
	echo ', "$(word 1, $^)" u 1:3 w lp ps 0.5 lc "blue"  dt 1 t "2nd" \$(void)' >> $(SCRIPT)
	echo ', "$(word 1, $^)" u 1:4 w lp ps 0.5 lc "black" dt 2 t "3rd" \$(void)' >> $(SCRIPT)
	echo ', "$(word 1, $^)" u 1:5 w lp ps 0.5 lc "red"   dt 1 t "4th" \$(void)' >> $(SCRIPT)
	echo ', "$(word 1, $^)" u 1:6 w lp ps 0.5 lc "black" dt 3 t "5th" \$(void)' >> $(SCRIPT)
	#echo ', "$(word 1, $^)" u 1:7 w lp lc "black" dt 6 t "6th" \$(void)' >> $(SCRIPT)
	$(PLOTTER) $(SCRIPT)
	gv $@
#==============================================================================#

task_a.txt energies.txt plot.txt: main
	./$<

main: main.c
#==============================================================================#

.PHONEY:clean
clean:                     # this is "clean" target. it does not depend on anything
	rm -f $(TARGETS)
	find ./ -type f -executable -delete
	$(RM) *.dat *.txt $(SCRIPT) *.svg *.png *.pdf log* *.eps


