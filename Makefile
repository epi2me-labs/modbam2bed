OS := $(shell uname)
ifeq ($(OS), Darwin)
SEDI=sed -i '.bak'
export LIBRARY_PATH=/usr/local/Cellar/openssl@1.1/1.1.1j/lib
else
SEDI=sed -i
endif


libhts.a:
	@echo Compiling $(@F)
	cd htslib/ && autoheader && autoconf && CFLAGS=-fpic ./configure && make
	cp htslib/$@ $@


.PHONY: clean_htslib
clean_htslib:
	cd htslib && make clean || exit 0


pileup: libhts.a src/medaka_common.c src/medaka_counts.c src/medaka_bamiter.c
	gcc -pthread  -g -Wall -fstack-protector-strong -D_FORTIFY_SOURCE=2 -fPIC -std=c99 -msse3 -O3 \
		-Isrc -Ihtslib \
		src/medaka_common.c src/medaka_counts.c src/medaka_bamiter.c libhts.a \
		-lm -lz -llzma -lbz2 -lpthread -lcurl -lcrypto \
		-o $(@)


.PHONY: clean
clean: clean_htslib
	rm pileup


