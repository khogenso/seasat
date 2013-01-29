target: all

INCLUDES = -Ilibsgp4 -I../include
LIBS = libfftw3f.a libsgp4.a
SRC = create_roi_in.c \
	dates.c \
	dop.c \
	spectra.c \
	fix_vectors.c \
	get_peg.c \
	runtest.cpp \
	tle_code.c
	
	
libsgp4.a:
	( cd libsgp4; make; mv libsgp4.a ..; cd .. )

all: libsgp4.a
	c++ -o create_roi_in $(SRC) $(INCLUDES) $(LIBS) -lm

clean:
	rm -f *.o create_roi_in
