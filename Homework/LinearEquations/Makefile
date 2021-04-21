cc = gcc
CFLAGS += -Wall -pedantic -std=gnu1x $$(gsl-config --cflags)
LDLIBS += $$(gsl-config --libs)

EPSABS = 1e-8
EPSREL = 1e-8
FILEPOINTER = output.txt

# For gnuplot
SCRIPT = log.gpi
FORMAT = eps
PLOTTER = gnuplot
void =


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
	for i in $$(seq 10 5 200);\
		do\
		time --format "$$i\t%R" --append --output $@ ./$< -n $$i -m $$i;\
		done

convergence_input.txt: my_timer
	for i in $$(seq 10 5 200);\
		do\
		time --format "$$i\t%R" --append --output $@ ./$< -n $$i -m $$i -a $(EPSABS) -r $(EPSREL) -f $(FILEPOINTER);\
		done


# Stripped down for timing only part a
gsl_timer: gsl_timer.c
my_timer: my_timer.c

prog: main main.c
	./$<
main: main.c


.PHONEY: clean
clean:
	find ./ -type f -executable -delete
	rm *.txt *.gpi *.log *.eps
