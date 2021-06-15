CC = gcc
LDLIBS = $$(gsl-config --libs)
CFLAGS = $$(gsl-config --cflags) -pedantic -Wall

OBJECTS = main contour.gpi contour.png dnewton.txt broyden.txt hybrid.txt hybrids.txt convergence.txt data2.txt task2.txt swave.png swave.gpi convergence.gpi convergence.png precision.txt

all: swave.gpi contour.gpi convergence.gpi

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




swave.gpi: main data2.txt Makefile
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
	echo 'plot "data2.txt" u 1:2 w l lw 2 dt 2 lc "blue" t "Estimate", \\' >> $@
	echo '"data2.txt" u 1:3      w l lw 2 dt 5 lc "red"   t "Actual"' >> $@
# Residual
	echo 'set tmargin 0.5; set bmargin 4' >> $@
	echo 'set xtics scale 0.5 font ".8" nomirror' >> $@
	echo 'set xlabel "r [Bohr Radii]"' >> $@
	echo 'plot "data2.txt" u 1:($$3-$$2) w l lt 1 t "Residual"' >> $@
# Multiplot done
	echo 'unset multiplot' >> $@
# Bash commands
	gnuplot $@
	ristretto swave.png

contour.gpi: main dnewton.txt broyden.txt hybrid.txt hybrids.txt Makefile
	echo 'set terminal pngcairo truecolor' > $@
	echo 'set output "contour.png"' >> $@
	echo 'set title "Contour of Rosenbrock Function"' >> $@
	echo 'set xlabel "x"' >> $@
	echo 'set ylabel "y"' >> $@
	echo 'set zlabel "z"' >> $@
	echo 'set xtics 1' >> $@
	echo 'set ytics 1' >> $@
	echo 'set key default' >> $@
	# View
	#echo 'set view 120, 210, 1, 1.2' >> $@
	echo 'set view 120, 30, 1, 1' >> $@
	# Contour plot
	echo 'set contour' >> $@
	echo 'set contour surface' >> $@
	echo 'set contour base' >> $@
	#echo 'set pm3d map' >> $@
	echo 'set isosamples 15, 15' >> $@
	# The color bar (left)
	echo 'set cntrparam level incremental 5e2, 5e3, 5e4' >> $@
#	echo 'set palette rgbformulae 33,13,10' >> $@
	# Margins
#	echo 'set lmargin at screen 0.1' >> $@
#	echo 'set bmargin at screen 0.1' >> $@
#	echo 'set rmargin at screen 0.80' >> $@
#	echo 'set tmargin at screen 0.85' >> $@
	# Plot functions
	echo 'Rosenbrock(x,y) = (1-x)**2 + 100*(y - x**2)**2' >> $@
	echo 'splot [-2:2] [-2:2.5] Rosenbrock(x,y) notitle ,\\' >> $@
	echo '"hybrids.txt" u "x":"y":"f(x,y)" w lp lw 1 dt 1 lc rgb "black"   t "hybrids",\\' >> $@
	echo '"broyden.txt" u "x":"y":"f(x,y)" w lp lw 1 dt 1 lc rgb "green"  t "broyden",\\' >> $@
	echo '"hybrid.txt"  u "x":"y":"f(x,y)" w lp lw 1 dt 1 lc rgb "blue" t "hybrid" ,\\' >> $@
	echo '"dnewton.txt" u "x":"y":"f(x,y)" w lp lw 1 dt 1 lc rgb "red" t "dnewton",\\' >> $@
	echo '"-" w p pt 5 ps 3 t "start"'>> $@
# We give in data by the following lines to the "-". e is the end of data.
	echo '-2.0 3.0 109.0' >> $@
	echo 'e' >> $@
	gnuplot $@
	ristretto contour.png

convergence.txt precision.txt hybrid.txt hybrids.txt dnewton.txt\
broyden.txt data2.txt: main
	./main

main: main.c


.PHONEY: clean
clean:
	rm $(OBJECTS)
