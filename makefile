flags=

ifeq ($(CC), clang)
ifeq ($(CONFIG), debug)
flags=-D_GNU_SOURCE -DDEBUG -O0 -g \
	-Weverything \
	-Wno-bad-function-cast -Wno-disabled-macro-expansion \
	-Wno-aggregate-return -Wno-covered-switch-default -Wno-unused-function
endif
ifeq ($(CONFIG), release)
flags=-DNDEBUG -Ofast -Wl,-s -w
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
flags=-D_GNU_SOURCE -DNDEBUG -Ofast -s -w
endif
endif

CFLAGS=-D_POSIX_C_SOURCE=200809L -std=c11 `pkg-config --cflags gsl` $(flags)
LDLIBS=-lm `pkg-config --libs gsl`

plot: plots.pdf

run: build
	GSL_RNG_TYPE=mt19937 GSL_RNG_SEED=0 time -v \
	./qho -d 1 -N 1 -M 32 -K 64 -h 8192 -p 131072 -H 16 -P 256 -T 0.1

check: build
	cppcheck -I/usr/include --enable=all *.c *.h
	valgrind --leak-check=full --tool=memcheck \
	./qho -d 2 -N 4 -M 8 -K 16 -h 128 -p 256 -H 32 -P 64 -T 0.1

build: he4 qho

clean: shallow-clean
	$(RM) he4 qho run-latest
	# TODO Remove the 2.
	$(RM) -r run-2*
	$(RM) *.eps *.tex plots.pdf

shallow-clean:
	$(RM) *.o
	$(RM) *-eps-converted-to.pdf *.aux *.log

plots.pdf: plots.tex
	pdflatex $<

plots.tex: energy.tex raddist.tex params.tex proj.tex \
	posdist-1.tex pots-1.tex \
	polys-2.tex posdist-2.tex pots-2.tex \
	polys-3.tex
	./plotgen $(basename $^) > $@

he4: he4.o err.o fp.o hist.o opt.o ran.o secs.o sigs.o sim.o size.o stats.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

qho: qho.o err.o fp.o hist.o opt.o ran.o secs.o sigs.o sim.o size.o stats.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

%.tex: %.gp
	gnuplot $<

%.o: %.c *.h
	$(CC) $(CFLAGS) -c -o $@ $<
