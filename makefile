flags=

ifeq ($(CC), clang)
ifeq ($(CONFIG), debug)
flags=-O0 -g -Weverything \
	-Wno-covered-switch-default -Wno-unused-macros
endif
ifeq ($(CONFIG), release)
flags=-Ofast -Wl,-s -w
endif
endif

ifeq ($(CC), gcc)
ifeq ($(CONFIG), debug)
flags=-O0 -g `cat gcc-$$(./gcc-version | tr . _)-release` \
	-Wno-error -Wno-fatal-errors -Wno-system-headers \
	-Wno-c++-compat -Wno-declaration-after-statement \
	-Wno-traditional -Wno-traditional-conversion \
	-Wno-unsuffixed-float-constants \
	-Wno-address
endif
ifeq ($(CONFIG), release)
flags=-Ofast -s -w
endif
endif

CFLAGS=-D_POSIX_C_SOURCE=200809L -std=c11 `pkg-config --cflags gsl` $(flags)
LDLIBS=-lm `pkg-config --libs gsl`

run: build
	./qho
	# gnuplot -p qho-ensemble-2d.gp
	gnuplot -p qho-ensemble-3d.gp

check: build
	valgrind --leak-check=full --tool=memcheck ./qho

build: qho

clean: shallow-clean
	$(RM) he4 qho

shallow-clean:
	$(RM) *.o

qho: floating.o nth.o qho.o report.o size.o timepack.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

%.o: %.c *.h
	$(CC) $(CFLAGS) -c -o $@ $<
