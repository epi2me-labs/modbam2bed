OS := $(shell uname)
ifeq ($(OS), Darwin)
    SEDI=sed -i '.bak'
    export LIBRARY_PATH=/usr/local/Cellar/openssl@1.1/1.1.1k/lib
    ARGP ?= /usr/local/Cellar/argp-standalone/1.3/lib/libargp.a
else
    SEDI=sed -i
    ARGP ?=
endif


libhts.a:
	@echo Compiling $(@F)
	cd htslib/ && autoheader && autoconf && CFLAGS=-fpic ./configure && make
	cp htslib/$@ $@


.PHONY: clean_htslib
clean_htslib:
	cd htslib && make clean || exit 0


modbam2bed: src/common.c src/counts.c src/bamiter.c src/args.c libhts.a
	gcc -pthread  -g -Wall -fstack-protector-strong -D_FORTIFY_SOURCE=2 -fPIC -std=c99 -msse3 -O3 \
		-Isrc -Ihtslib \
		$^ $(ARGP) \
		-lm -lz -llzma -lbz2 -lpthread -lcurl -lcrypto \
		-o $(@)


.PHONY: clean
clean: clean_htslib
	rm modbam2bed


