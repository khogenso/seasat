target: all

all:
	make -C create_roi_in
	make -C fix_headers
	make -C decoder
	cp create_roi_in/create_roi_in bin
	cp create_roi_in/SEASAT_TLEs.txt bin
	cp fix_headers/fix_headers bin
	cp fix_headers/fix_stairs bin
	cp fix_headers/dis_search bin
	cp decoder/seasat_decoder bin

clean:
	make clean -C create_roi_in
	make clean -C fix_headers
	make clean -C decoder
	rm -f bin/*
