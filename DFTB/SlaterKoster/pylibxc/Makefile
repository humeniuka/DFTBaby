# path to the headers and library file of libxc 
# are assumed to be set in environment variables CPLUS_INCLUDE_PATH and LD_LIBRARY_PATH
# libxc must have been compiled as a shared library using
#  ./configure --enable-shared
all: _pylibxc.so

# ---- Link ---------------------------
_pylibxc.so:  pylibxc.o
	g++ -Wall -shared -o _pylibxc.so pylibxc.o -lxc


# ---- gcc C compile ------------------
INC=/usr/include/python2.7
pylibxc.o:  pylibxc.cpp
	g++ -Wall -fPIC -c pylibxc.cpp -I $(INC) -I $(INC)/numpy

# Clean
clean:
	rm pylibxc.o
	rm _pylibxc.so

