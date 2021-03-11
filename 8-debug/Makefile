CFLAGS = -o -w -std=gnu99
CFLAGS += $$(gsl-config --cflags) # Changed
LDLIBS += $$(gsl-config --libs)

all: output ; cat $< # Changed $@ (which refers to all) to $< (output)
output: main ; ./$< > $@ # Changed $@ to ./$< and $< to $@
clean: ; $(RM) main output #Changed $RM to $(RM)

