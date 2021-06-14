# Compiler
CC = gcc
CFLAGS += -Wall -pedantic -fwrapv -O3 -std=gnu1x $$(gsl-config --cflags)
LDLIBS += $$(gsl-config --libs)
TARGETS = results.txt log.gpi operations.svg convergence_input.txt convergence_gsl_input.txt gsl.txt

#==============================================================================#
# For gnuplot
SCRIPT = log.gpi
FORMAT = svg
PLOTTER = gnuplot
void =
#==============================================================================#
# CONVERGENCE TESTS:
# Figure
operations.$(FORMAT): convergence_input.txt convergence_gsl_input.txt Makefile
	echo 'set terminal eps' > $(SCRIPT)
	echo 'set output "$@"' >> $(SCRIPT)
	echo 'set key left' >> $(SCRIPT)
	echo 'set title "Convergence of SVG alg. of square matrices"' >> $(SCRIPT)
	echo 'set grid' >> $(SCRIPT)
	echo 'set xlabel "n [Dimension of square matrix]"' >> $(SCRIPT)
	echo 'set ylabel "t [Some time unit]"' >> $(SCRIPT)
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

convergence_input.txt: main
	for i in $$(seq 20 5 200);\
		do\
		time --format "$$i\t%R" --append --output $@ ./$< -n $$i -m $$i;\
		done

gsl.txt: gsl_timer
	./$<

gsl_timer: gsl_timer.c


results.txt: main
	./$<

main: main.c

#==============================================================================#

.PHONEY:clean
clean:                     # this is "clean" target. it does not depend on anything
	find ./ -type f -executable -delete
	rm $(TARGETS)


