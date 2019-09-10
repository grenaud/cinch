
all: lib/insertSize/src/insize

lib/insertSize/src/insize:
	mkdir lib/insertsize/
	git clone --depth 1 https://github.com/grenaud/insertsize.git lib/insertsize/
	make -C  lib/insertsize/src/

clean :
	make -C lib/insertsize/src/ clean

