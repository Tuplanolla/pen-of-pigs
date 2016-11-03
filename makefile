configflags=

ifeq ($(CC), clang)
ifeq ($(CONFIG), debug)
configflags=-O0 -g -Weverything \
	-Wno-covered-switch-default -Wno-unused-macros
endif
ifeq ($(CONFIG), release)
configflags=-Ofast -Wl,-s -w
endif
endif

ifeq ($(CC), gcc)
ifeq ($(CONFIG), debug)
configflags=-O0 -g `cat gcc-$$(./gcc-version | tr . _)-release` \
	-Wno-error -Wno-fatal-errors -Wno-system-headers \
	-Wno-c++-compat -Wno-declaration-after-statement \
	-Wno-traditional -Wno-traditional-conversion \
	-Wno-unsuffixed-float-constants \
	-Wno-address
endif
ifeq ($(CONFIG), release)
configflags=-Ofast -s -w
endif
endif

jobflags=

ifeq ($(JOB), qho)
jobflags=-DPEN_OF_QHO
endif

ifeq ($(JOB), he4)
jobflags=-DPEN_OF_HE4
endif

CFLAGS=-D_POSIX_C_SOURCE=200809L -std=c11 `pkg-config --cflags gsl` \
	$(configflags) $(jobflags)
LDLIBS=-lm `pkg-config --libs gsl`

run: build
	./pen-of-pigs

check: build
	valgrind --leak-check=full --tool=memcheck ./pen-of-pigs

build: pen-of-pigs

clean: shallow-clean
	$(RM) pen-of-pigs

shallow-clean:
	$(RM) *.o

pen-of-pigs: main.o floating.o report.o size.o timepack.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

%.o: %.c *.h
	$(CC) $(CFLAGS) -c -o $@ $<
