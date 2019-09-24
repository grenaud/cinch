all: lib/insertsize/src/insize  lib/hmmcopy_utils/bin/readCounter

lib/insertsize/src/insize.cpp:
	rm -rfv lib/insertsize/ && mkdir -p lib/insertsize/
	git clone --depth 1 https://github.com/grenaud/insertsize.git lib/insertsize/

lib/insertsize/src/insize: lib/insertsize/src/insize.cpp
	make -C  lib/insertsize/src/

lib/hmmcopy_utils/src/bin/readCounter.cpp:
	rm -rfv lib/hmmcopy_utils/ && mkdir -p lib/hmmcopy_utils/
	git clone --depth 1 https://github.com/shahcompbio/hmmcopy_utils.git  lib/hmmcopy_utils/

lib/hmmcopy_utils/bin/readCounter: lib/hmmcopy_utils/src/bin/readCounter.cpp
	cd  lib/hmmcopy_utils/ && cmake . && make && cd ../..

clean :
	make -C lib/insertsize/src/ clean
	make -C lib/hmmcopy_utils/ clean

