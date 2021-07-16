# Builds a cache of binaries which can just be copied for CI
OS := $(shell uname)
ifeq ($(OS), Darwin)
SEDI=sed -i '.bak'
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
		-o $(@) -std=c99 -msse3 -O3


modpileup: libhts.a pileup_mod.c
	gcc -pthread  -g -Wall -fstack-protector-strong -D_FORTIFY_SOURCE=2 -fPIC -std=c99 -msse3 -O3 \
		-Isrc -Ihtslib \
		pileup_mod.c libhts.a \
		-lm -lz -llzma -lbz2 -lpthread -lcurl -lcrypto \
		-o $(@) -std=c99 -msse3 -O3
