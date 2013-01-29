target: all

all: fix_headers fix_stairs dis_search

fix_headers:
	c++ -o fix_headers fix_headers.c -lm

fix_stairs:
	c++ -o fix_stairs fix_stairs.c -lm

dis_search:
	c++ -o dis_search search.c -lm

clean:
	rm -f dis_search fix_stairs fix_headers

