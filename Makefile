OS := $(shell uname)
ifeq ($(OS), Darwin)
    SEDI=sed -i '.bak'
	# mainly for dev builds using homebrew things
    export LIBRARY_PATH=/usr/local/Cellar/openssl@1.1/1.1.1k/lib
    ARGP ?= /usr/local/Cellar/argp-standalone/1.3/lib/libargp.a
else
    SEDI=sed -i
    ARGP ?=
endif

CC ?= gcc
CFLAGS ?= -fpic -msse3 -O3
EXTRA_CFLAGS ?=
EXTRA_LDFLAGS ?=
EXTRA_LIBS ?=
HTS_CONF_ARGS ?=


.PHONY: default
default: modbam2bed


htslib/libhts.a:
	@echo Compiling $(@F)
	cd htslib/ \
		&& autoheader \
		&& autoconf \
		&& CFLAGS="$(CFLAGS)" ./configure $(HTS_CONF_ARGS) \
		&& make -j 4


.PHONY: clean_htslib
clean_htslib:
	cd htslib && make clean || exit 0


modbam2bed: src/common.c src/counts.c src/bamiter.c src/args.c htslib/libhts.a
	$(CC) -pthread -Wall -fstack-protector-strong -D_FORTIFY_SOURCE=2 $(CFLAGS) \
		-Isrc -Ihtslib $(EXTRA_CFLAGS) $(EXTRA_LDFLAGS)\
		$^ $(ARGP) \
		-lm -lz -llzma -lbz2 -lpthread -lcurl -lcrypto $(EXTRA_LIBS) \
		-o $(@)


.PHONY: clean
clean: clean_htslib
	rm -rf modbam2bed


.PHONY: mem_check
mem_check: modbam2bed
	valgrind --error-exitcode=1 --tool=memcheck --leak-check=full --show-leak-kinds=all -s \
		./modbam2bed -b 0.66 -a 0.33 -t 2 -r ecoli1 test_data/400ecoli.bam test_data/ecoli.fasta > /dev/null


