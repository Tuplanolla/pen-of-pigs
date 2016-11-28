flags=

ifeq ($(CC), clang)
ifeq ($(CONFIG), debug)
flags=-DDEBUG -O0 -g -Weverything \
	-Wno-covered-switch-default -Wno-unused-macros
endif
ifeq ($(CONFIG), release)
flags=-DNDEBUG -Ofast -Wl,-s -w
endif
endif

ifeq ($(CC), gcc)
ifeq ($(CONFIG), debug)
flags=-DDEBUG -O0 -g `cat gcc-$$(./gcc-version | tr . _)-release` \
	-Wno-error -Wno-fatal-errors -Wno-system-headers \
	-Wno-c++-compat -Wno-declaration-after-statement \
	-Wno-traditional -Wno-traditional-conversion \
	-Wno-unsuffixed-float-constants \
	-Wno-address
endif
ifeq ($(CONFIG), release)
flags=-DNDEBUG -Ofast -s -w
endif
endif

flags+=-Wno-unused-function # TODO Remove this.

CFLAGS=-D_POSIX_C_SOURCE=200809L -std=c11 `pkg-config --cflags gsl` $(flags)
LDLIBS=-lm `pkg-config --libs gsl`

run: build
	GSL_RNG_TYPE=mt19937 GSL_RNG_SEED=0 time -v ./qho

check: build
	valgrind --leak-check=full --tool=memcheck ./qho

build: qho

clean: shallow-clean
	$(RM) he4 qho

shallow-clean:
	$(RM) *.o

qho: err.o fp.o nth.o phys.o qho.o secs.o sigs.o size.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

%.o: %.c *.h
	$(CC) $(CFLAGS) -c -o $@ $<
