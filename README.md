# Pen of PIGS

You can compile and run the program as follows.

    $ make clean build
    $ make run

Options exist.

    $ make CC=clang TARGET=debug JOB=qho clean build

The available targets are

* `run` to run the program,
* `check` to check the program for memory leaks,
* `build` to compile the program,
* `clean` to remove all generated files and
* `shallow-clean` to remove only unimportant generated files.

The supported compilers (`CC`) are

* `clang` and
* `gcc`.

The possible configurations (`CONFIG`) are

* `debug` and
* `release`.
