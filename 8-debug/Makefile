CFLAGS = -o -w -std=gnu99

# Changed FLAGS and LIBS from single quotation (weak), to correct Make syntax
CFLAGS += $$(gsl-config --cflags) # Changed
LDLIBS += $$(gsl-config --libs)

# Changed $@ (which refers to all) to $< (output)
all: output ; cat $<


# Changed $@ to ./$< and $< to $@
output: main ; ./$< > $@

#Changed $RM to $(RM), since $ would act on R only if no parenthesis
clean: ; $(RM) main output

