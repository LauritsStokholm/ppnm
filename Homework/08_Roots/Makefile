# Compiler
CC = gcc
CFLAGS += -Wall -pedantic -std=gnu1x $$(gsl-config --cflags)
LDLIBS += $$(gsl-config --libs)

#==============================================================================#
# For gnuplot
SCRIPT = log.gpi
FORMAT = png
PLOTTER = gnuplot
void =


tagets: convergence.gpi swave.gpi
#==============================================================================#
convergence.gpi: main convergence.txt precision.txt Makefile
	echo 'set terminal pngcairo truecolor' > $@
	echo 'set output "convergence.png"' >> $@
	echo 'set title "Exercise b and c) Convergence in Shooting Method"' >> $@
	echo 'set key bottom right' >> $@
	echo 'set xlabel "Value of rmax"' >>$@
	echo 'set ylabel "Number of Itterations for convergence"' >>$@
	echo 'plot "convergence.txt" u 1:2 w l lw 2 dt 2 lc "blue" t "Boundary Condition f(rmax)=0",\\' >> $@
	echo '"precision.txt" u 1:2 w l lw 2 dt 4 lc "red" t "Boundary Condition f(rmax)=r*exp(-kr)"' >> $@
	gnuplot $@
	ristretto convergence.png


# Figure 1
swave.gpi: task_b.dat Makefile
	echo 'set terminal pngcairo truecolor' > $@
	echo 'set output "swave.png"' >> $@
# Multiplot settings
	echo 'set tmargin 1; set bmargin 0.5' >> $@
	echo 'set lmargin 8; set rmargin 5' >> $@
	echo 'unset xtics' >> $@
	echo 'set ytics nomirror' >> $@
	echo 'set ylabel "psi(r)"' >> $@
	echo 'set multiplot layout 2,1 title "Solution to s-wave\\n"' >> $@
# Data plot
	echo 'plot "$<" u 1:2 w l lw 2 dt 2 lc "blue" t "Estimate", \\' >> $@
	echo '     "$<" u 1:3 w l lw 2 dt 5 lc "red"  t "Actual"' >> $@
# Residual
	echo 'set tmargin 0.5; set bmargin 4' >> $@
	echo 'set xtics scale 0.5 font ".8" nomirror' >> $@
	echo 'set xlabel "r [Bohr Radii]"' >> $@
	echo 'plot "$<" u 1:($$3-$$2) w l lt 1 t "Residual"' >> $@
# Multiplot done
	echo 'unset multiplot' >> $@
# Bash commands
	gnuplot $@
	ristretto swave.png


#	echo 'set terminal svg' > $(SCRIPT)
#	echo 'set output "$@"' >> $(SCRIPT)
#	echo 'set key out' >> $(SCRIPT)
#	echo 'set title "Solution to s-wave for H-atom"' >> $(SCRIPT)
#	echo 'set grid' >> $(SCRIPT)
#	echo 'set xlabel "r [Bohr radii]"' >> $(SCRIPT)
#	echo 'set ylabel "psi(r)"' >> $(SCRIPT)
#	echo 'plot \$(void)' >> $(SCRIPT)
#	echo ' "$(word 1, $^)" u 1:2 w l dt 1 lc "blue" t "numerical" \$(void)' >> $(SCRIPT)
#	echo ',"$(word 1, $^)" u 1:3 w l dt 2 lc "red" t "expected" \$(void)' >> $(SCRIPT)
##	echo ' "$(word 1, $^)" u 1:2 w l dt 1 lc "black" t "rk12, y0" \$(void)' >> $(SCRIPT)
##	echo ',"$(word 1, $^)" u 1:3 w l dt 1 lc "black" t "rk12, y1" \$(void)' >> $(SCRIPT)
##	echo ',"$(word 1, $^)" u 1:4 w l dt 2 lc "red"   t "cos" \$(void)' >> $(SCRIPT)
##	echo ',"$(word 1, $^)" u 1:5 w l dt 2 lc "blue"  t "sin" \$(void)' >> $(SCRIPT)
#	$(PLOTTER) $(SCRIPT)
#	ristretto $@
#==============================================================================#

convergence.txt precision.txt task_b.dat task_c.dat: main
	./main

main: main.c


.PHONEY:clean
clean:                     # this is "clean" target. it does not depend on anything
	rm -f $(TARGETS)
	find ./ -type f -executable -delete
	$(RM) *.txt *.dat $(SCRIPT) *.svg *.png *.pdf log* *.eps *.gpi


