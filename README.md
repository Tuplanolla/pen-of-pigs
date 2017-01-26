# Pen of PIGS

This repository contains a simplified implementation
of the path integral ground state quantum Monte Carlo method.
While the code seems to produce correct results and run reasonably fast,
it is primarily intended for pedagogical purposes.
Clarity and conciseness of the presentation are favored
over numerical stability and performance tricks.
If you would prefer something with the opposite design goals,
see [WA-PIMC][wa-pimc] by Adrian Del Maestro instead.

## Overview

Pen of PIGS was written by Sampsa "Tuplanolla" Kiiskinen
to support his bachelor's thesis on quantum statistical mechanics.
Even though the project was in the works for a while,
the bulk of the work happened between 2016-11-01 and 2017-01-23.

## License

Pen of PIGS is free software and as such
licensed under the GNU General Public License version 3 or later.
The full license can be found in the `LICENSE` file that
resides in the same directory as this file.
In short, copies and derivative works are permitted
as long as they use a compatible license.

## Requirements

Pen of PIGS is written in C (per ISO/IEC 9899:2011) and
makes heavy use of POSIX features (per POSIX.1-2008).
Its only dependency, aside from build essentials,
is the [GNU Scientific Library][gsl]
(available, at least in the case of [Debian][debian],
in the `libgsl0-dev` package).

Although graphics support is not necessary,
a handful of useful visualizations are provided as [Gnuplot][gnuplot] scripts,
so it may be a good idea to install that as well
(the scripts use the `wxt` and `epslatex` terminals by default,
which, for Debian, requires the `gnuplot-x11` package).

Debugging is not necessary either, but, if desired,
requires [Valgrind][valgrind] and [Cppcheck][cppcheck].

Some automation is also available in the form of [Bourne shell][sh] scripts.

## Usage

Pen of PIGS comes with several frontends for different kinds of computations.
You can compile and run them as follows.

    $ make build
    $ make run

Note that these invocations use the default build environment,
which depends on your system and, while reliable, may not be ideal.
For this reason specialized options and targets exist.

The available targets are

* `run` to run the program,
* `check` to check the program for memory leaks and other potential problems,
* `build` to compile the program,
* `plot` to render visualizations,
* `clean` to remove all generated files and
* `shallow-clean` to remove only unimportant generated files.

The supported compilers (`CC`) are

* `clang` for the [C Language Family Frontend for LLVM][clang] and
* `gcc` for the [GNU Compiler Collection][gcc].

The possible configurations (`CONFIG`) are

* `debug` for debugging and
* `release` for real computations.

You can put these into good use as follows.

    $ make CC=clang CONFIG=debug clean check
    $ make CC=gcc CONFIG=release clean run

## Anything Else

Read the source code.
It is not scary.
I promise.

[wa-pimc]: http://code.delmaestro.org/
[gsl]: https://www.gnu.org/software/gsl/
[debian]: https://www.debian.org/
[gnuplot]: http://gnuplot.sourceforge.net/
[valgrind]: http://valgrind.org/
[cppcheck]: http://cppcheck.sourceforge.net/
[sh]: http://gondor.apana.org.au/~herbert/dash/
[clang]: http://clang.llvm.org/
[gcc]: https://gcc.gnu.org/
