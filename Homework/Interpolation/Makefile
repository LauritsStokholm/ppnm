CC = gcc
CFLAGS += -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS += $$(gsl-config --libs)

# For gnuplot
SCRIPT = log.gpi
FORMAT = eps
PLOTTER = gnuplot
void =

TARGETS = data.txt main lspline.$(FORMAT) cspline.$(FORMAT)


default: $(TARGETS)

# Figure 2
cspline.$(FORMAT): data.txt cspline.txt cgsl.txt Makefile
	echo 'set terminal eps' > $(SCRIPT)
	echo 'set output "$@"' >> $(SCRIPT)
	echo 'set key left' >> $(SCRIPT)
	echo 'set title "This is a title"' >> $(SCRIPT)
	echo 'set grid' >> $(SCRIPT)
	echo 'set xlabel "x"' >> $(SCRIPT)
	echo 'set ylabel "y"' >> $(SCRIPT)
	echo 'plot \$(void)' >> $(SCRIPT)
	echo '"$(word 1, $^)" u 1:2 w p pt 2 ps 1 t "Data"  \$(void)' >> $(SCRIPT)
	echo ', "$(word 2, $^)" u 1:2 w l dt 1 lc "black" lw 5 t "cspline_eval"\$(void)' >> $(SCRIPT)
	echo ', "$(word 2, $^)" u 1:3 w l dt 1 lc "red" lw 5 t "cspline_deriv"\$(void)' >> $(SCRIPT)
	echo ', "$(word 2, $^)" u 1:4 w l dt 1 lc "red" lw 5 t "cspline_integ"\$(void)' >> $(SCRIPT)
	echo ', "$(word 3, $^)" u 1:2 w l dt 2 lc "blue"  lw 5 t "gsl cspline_eval"\$(void)' >> $(SCRIPT)
	echo ', "$(word 3, $^)" u 1:3 w l dt 2 lc "green"  lw 5 t "gsl cspline_deriv"\$(void)' >> $(SCRIPT)
	echo ', "$(word 3, $^)" u 1:4 w l dt 2 lc "grey"  lw 5 t "gsl cspline_integ"\$(void)' >> $(SCRIPT)
	$(PLOTTER) $(SCRIPT)
	gv $@

# Figure 1
lspline.$(FORMAT): data.txt lspline.txt lgsl.txt qspline.txt Makefile
	echo 'set terminal eps' > $(SCRIPT)
	echo 'set output "$@"' >> $(SCRIPT)
	echo 'set key left' >> $(SCRIPT)
	echo 'set title "This is a title"' >> $(SCRIPT)
	echo 'set grid' >> $(SCRIPT)
	echo 'set xlabel "x"' >> $(SCRIPT)
	echo 'set ylabel "y"' >> $(SCRIPT)
	echo 'plot \$(void)' >> $(SCRIPT)
	echo '"$(word 1, $^)" u 1:2 w p pt 2 ps 1 t "Data"  \$(void)' >> $(SCRIPT)
	echo ', "$(word 2, $^)" u 1:2 w l dt 1 lc "black" lw 5 t "linterp_eval"\$(void)' >> $(SCRIPT)
	echo ', "$(word 2, $^)" u 1:3 w l dt 1 lc "black" lw 5 t "linterp_deriv"\$(void)' >> $(SCRIPT)
	echo ', "$(word 2, $^)" u 1:4 w l dt 1 lc "black" lw 5 t "linterp_integ"\$(void)' >> $(SCRIPT)
#
	echo ', "$(word 3, $^)" u 1:2 w l dt 2 lc "blue"  lw 5 t "gsl linterp_eval"\$(void)' >> $(SCRIPT)
	echo ', "$(word 3, $^)" u 1:3 w l dt 2 lc "blue"  lw 5 t "gsl linterp_deriv"\$(void)' >> $(SCRIPT)
	echo ', "$(word 3, $^)" u 1:4 w l dt 2 lc "blue"  lw 5 t "gsl linterp_integ"\$(void)' >> $(SCRIPT)
#
	echo ', "$(word 4, $^)" u 1:2 w l dt 3 lc "yellow"  lw 5 t "quadratic eval"\$(void)' >> $(SCRIPT)
	echo ', "$(word 4, $^)" u 1:3 w l dt 3 lc "red"  lw 5 t "quadratic deriv"\$(void)' >> $(SCRIPT)
	echo ', "$(word 4, $^)" u 1:4 w l dt 3 lc "yellow"  lw 5 t "quadratic integ"\$(void)' >> $(SCRIPT)
	$(PLOTTER) $(SCRIPT)
	gv $@


lspline.txt lgsl.txt: main data.txt
	./$<

main: main.c

data.txt:
	for x in `seq 0 1 10`; do echo "$$x\t$$(($$x * $$x * $$x + $$x * $$x - $$x - 2))" >> $@; done;



.PHONEY:clean
clean:                     # this is "clean" target. it does not depend on anything
	find ./ -type f -executable -delete
	$(RM) *.txt *.dat $(SCRIPT) *.svg *.png *.pdf log* *.eps


