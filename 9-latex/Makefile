SHELL = /bin/bash # Ubuntu uses dash, which seem to make trouble with interpretations of double backslash (see https://github.com/LauritsStokholm/makeitwork)
CC = gcc
CFLAGS = -Wall
LDLIBS = -lm

BLANK :=
SLASH = \$(BLANK)
HEADER = $(SLASH)$(SLASH)bfseries #and other settings
OBJECTS = exp.tex output.txt
LATEXTERMINAL = cairolatex pdf input blacktext header "$(HEADER)" font ",12" size 3.06, 3.06


main.pdf: main.tex exp.tex
	pdflatex $<
	gv $@

# To know size (given in inches pr. default/ can be set to cm)
# I have put \showthe\columnwidth at place in .tex file for plot to be,
# which stops compilation and prints columnwidth in pts. This was converted to
# inches, and used here

exp.tex: output.txt Makefile
	echo '\
	set terminal $(LATEXTERMINAL);\
	set output "$@";\
	set xlabel "x"; set ylabel "y";\
	f(x) = exp(x);\
	set tmargin 1; set bmargin 0.5;\
	set lmargin 8; set rmargin 5;\
	unset xtics;\
	set ytics nomirror;\
	set multiplot layout 2,1 title "Implementation of exponential function";\
	plot \
	  "$<" u "x":"myexp" w l lt 1 lw 7 dt 1 lc "black" t "self-implemented"\
	, "$<" u "x":"exp"   w l lt 2 lw 7 dt 2 lc "red" t "math.h"\
	, f(x) w l lt 2 lw 7 dt 3 lc "yellow" t "gnuplot";\
	set tmargin 0.5; set bmargin 4;\
	set yrange [-1e-11:1e-11];\
	set xtics scale 0.5 font ".8" nomirror;\
	set xlabel "x";\
	plot \
	"$<" u 1:($$3-$$2) w l lt 1 t "Residual";\
	' | gnuplot

output.txt: main
	./main

main: main.c

clean:
	$(RM) *.aux *.log *.pdf *.fls *.dvi *.fdb_latexmk *.synctex.gz *.bcf *.eps *.out *.run.xml *.pyg *.bbl *.blg
	$(RM) $(OBJECTS)
	find ./ -type f -executable -delete

