COMPILER = gcc
FLAGS = -O2 -Wall
LIBS = -lm -llapack

OBJECTS = ff.o linked_list.o input.o

../ff.so: $(OBJECTS) ff_pythonmodule.c
	rm -rf build
	python setup.py build
	mv build/lib.linux-x86_64-2.7/ff.so $@
	rm -rf build

clean:
	rm -rf build
	rm -f *.o

linked_list.o: linked_list.h
input.o: input.h atoms.h dreiding.h
ff.o: ff.h
