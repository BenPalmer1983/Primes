.PHONY: all run clean

ALL=prime-bitarray prime-8bit prime-logical-array prime-bitarray-parallel

FORTRAN=gfortran
FORTRANFLAGS=-Ofast -march=native

all: $(ALL)

%: %.f08
	$(FORTRAN) $(FORTRANFLAGS) -o $@ $<

prime-bitarray-parallel: prime-bitarray-parallel.f08
	$(FORTRAN) $(FORTRANFLAGS) -fopenmp -o $@ $<

run: all
	@for program in $(ALL); do ./$${program}; done

clean:
	rm -f $(ALL) *.mod
