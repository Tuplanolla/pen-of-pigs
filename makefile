flags=

ifeq ($(CC), clang)
ifeq ($(CONFIG), debug)
flags=-D_GNU_SOURCE -DDEBUG -O0 -g \
	-Weverything \
	-Wno-bad-function-cast -Wno-disabled-macro-expansion \
	-Wno-aggregate-return -Wno-covered-switch-default -Wno-unused-function
endif
ifeq ($(CONFIG), release)
flags=-DNDEBUG -O3 -Wl,-s -w
endif
endif

ifeq ($(CC), gcc)
ifeq ($(CONFIG), debug)
flags=-D_GNU_SOURCE -DDEBUG -Og -g \
	`cat gcc-$$(./gcc-version | tr . _)-release` \
	-Wno-error -Wno-fatal-errors -Wno-system-headers \
	-Wno-c++-compat -Wno-declaration-after-statement \
	-Wno-traditional -Wno-traditional-conversion \
	-Wno-unsuffixed-float-constants \
	-Wno-address -Wno-bad-function-cast -Wno-long-long \
	-Wno-aggregate-return -Wno-switch-default -Wno-unused-function
endif
ifeq ($(CONFIG), release)
flags=-D_GNU_SOURCE -DNDEBUG -O3 -s -w
endif
endif

CFLAGS=-D_POSIX_C_SOURCE=200809L -std=c11 `pkg-config --cflags gsl` $(flags)
LDLIBS=-lm -lrt `pkg-config --libs gsl`

build: pigs pimc

plot: plots.pdf

run: build
	GSL_RNG_TYPE=mt19937 GSL_RNG_SEED=0 time -v ./pimc \
	-s qho -d 1 -N 1 -M 32 -k 64 -K 64 -h 8192 -p 131072 -H 16 -P 256 \
	-L 8.0 -m 1.0 -T 0.125

check: build
	cppcheck -I/usr/include --enable=all *.c *.h
	valgrind --leak-check=full --tool=memcheck ./pimc \
	-s qho -d 2 -N 8 -M 16 -k 4 -K 32 -h 256 -p 512 -H 64 -P 128 \
	-L 8.0 -m 1.0 -T 0.125

clean: shallow-clean
	$(RM) pigs pimc run-latest
	$(RM) -r run-*
	$(RM) *.eps *.tex plots.pdf

shallow-clean:
	$(RM) *.o
	$(RM) *-eps-converted-to.pdf *.aux *.log

plots.pdf: plots.tex
	pdflatex $<

plots.tex: energy.tex energy-bead.tex raddist.tex params.tex proj.tex \
	posdist-1.tex pots-1.tex \
	polys-2.tex posdist-2.tex pots-2.tex \
	polys-3.tex
	./plotgen $(basename $^) > $@

pigs: pigs.o err.o fp.o hist.o opt.o ran.o secs.o sigs.o sim.o size.o stats.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

pimc: pimc.o err.o fp.o hist.o opt.o ran.o secs.o sigs.o sim.o size.o stats.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

%.tex: %.gp
	gnuplot $<

%.o: %.c *.h
	$(CC) $(CFLAGS) -c -o $@ $<
