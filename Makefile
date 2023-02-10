OS := $(shell uname)
ARCH := $(shell arch)

OS := $(shell uname)
ifeq ($(OS), Darwin)
    # mainly for dev builds using homebrew things
    EXTRA_LDFLAGS ?= -L$(shell brew --prefix openssl@1.1)/lib
    ARGP ?= $(shell brew --prefix argp-standalone)/lib/libargp.a
    ARGP_INCLUDE ?= -I$(shell brew --prefix argp-standalone)/include
else
    ARGP ?=
    ARGP_INCLUDE ?=
endif


CC ?= gcc
CFLAGS ?= -fpic -msse3 -O3 -std=c99
DEFLATE ?= $(PWD)/libdeflate
STATIC_HTSLIB ?= htslib/libhts.a
EXTRA_CFLAGS ?=
EXTRA_LDFLAGS ?=
EXTRA_LIBS ?=
HTS_CONF_ARGS ?=
HTS_CONF_ENV ?= CFLAGS="$(CFLAGS) $(EXTRA_CFLAGS)"

WITHDEFLATE ?= 
DEFLATEREQ =
ifeq ($(WITHDEFLATE), 1)
CFLAGS += -I$(DEFLATE) -L$(DEFLATE)
HTS_CONF_ARGS += --with-libdeflate
HTS_CONF_ENV += LDFLAGS="-L$(DEFLATE)"
EXTRA_LIBS += -ldeflate
DEFLATEREQ = libdeflate/libdeflate.so.0
endif

NOTHREADS ?=
ifeq ($(NOTHREADS), 1)
	CFLAGS += -DNOTHREADS
endif

VALGRIND ?= valgrind


.PHONY: default
default: modbam2bed

libdeflate/libdeflate.so.0:
	@echo Compiling $(@F)
	cd libdeflate && make
	

htslib/libhts.a: $(DEFLATEREQ)
	@echo Compiling $(@F)
	cd htslib/ \
		&& autoreconf -i \
		&& autoheader \
		&& autoconf \
		&& $(HTS_CONF_ENV) ./configure $(HTS_CONF_ARGS) \
		&& make -j 4


.PHONY: clean_htslib
clean_htslib:
	rm -rf htslib/autom4te.cache/ 
	cd htslib && make clean || exit 0


%.o: src/%.c
	mkdir -p obj && \
		$(CC) -c -pthread -Wall -fstack-protector-strong -D_FORTIFY_SOURCE=2 $(CFLAGS) \
		-Isrc -Ihtslib $(ARGP_INCLUDE) $(EXTRA_CFLAGS) $^ -o $@

.PHONY: clean_obj
clean_obj:
	rm -rf *.o


modbam2bed: modbam2bed.o common.o counts.o bamiter.o args.o $(STATIC_HTSLIB)
	$(CC) -pthread -Wall -fstack-protector-strong -D_FORTIFY_SOURCE=2 $(CFLAGS) \
		-Isrc -Ihtslib $(EXTRA_CFLAGS) $(EXTRA_LDFLAGS)\
		$^ $(ARGP) \
		-lm -lz -llzma -lbz2 -lpthread -lcurl -lcrypto $(EXTRA_LIBS) \
		-o $(@)

.PHONY: clean
clean: clean_obj clean_htslib
	rm -rf modbam2bed

.PHONY: mem_check
mem_check: modbam2bed
	$(VALGRIND) --error-exitcode=1 --tool=memcheck --leak-check=full --show-leak-kinds=all -s \
		./modbam2bed --threshold 0.66 -t 2 -r ecoli1 test_data/ecoli.fasta.gz test_data/400ecoli.bam test_data/400ecoli.bam > /dev/null


.PHONY: test_api
test_api: python
	${IN_VENV} && pip install pytest
	${IN_VENV} && pytest test --doctest-modules

### Python

PYTHON ?= python3
VENV ?= venv
venv: ${VENV}/bin/activate
IN_VENV=. ./${VENV}/bin/activate

$(VENV)/bin/activate:
	test -d $(VENV) || $(PYTHON) -m venv $(VENV) --prompt "modbam"
	${IN_VENV} && pip install pip --upgrade
	${IN_VENV} && pip install setuptools

.PHONY: python
python: htslib/libhts.a pymod.a $(VENV)/bin/activate
	${IN_VENV} && pip install -r requirements.txt
	${IN_VENV} && WITHDEFLATE=$(WITHDEFLATE) LDFLAGS=$(EXTRA_LDFLAGS) pip install -e .

.PHONY: clean_python
clean_python: clean_obj
	rm -rf dist build modbampy.egg-info pymod.a libmodbampy.abi3.so ${VENV}

pymod.a: common.o bamiter.o counts.o
	ar rcs $@ $^

test_python: python
	${IN_VENV} && pip install flake8 flake8-rst-docstrings flake8-docstrings flake8-import-order
	${IN_VENV} && flake8 modbampy \
		--import-order-style google --application-import-names modbampy,libmodbampy \
		--statistics
	${IN_VENV} && modbampy test_data/400ecoli.bam ecoli1 0 4000000 | wc -l
	${IN_VENV} && modbampy test_data/400ecoli.bam ecoli1 0 4000000 --pileup | wc -l

IN_BUILD=. ./pypi_build/bin/activate
pypi_build/bin/activate:
	test -d pypi_build || $(PYTHON) -m venv pypi_build --prompt "(pypi) "
	${IN_BUILD} && pip install pip --upgrade
	${IN_BUILD} && pip install --upgrade pip setuptools twine wheel readme_renderer[md] keyrings.alt

.PHONY: sdist
sdist: pypi_build/bin/activate
	${IN_BUILD} && python setup.py sdist


.PHONY: wheels
wheels: clean clean_python 
	docker run -v `pwd`:/io quay.io/pypa/manylinux2010_x86_64 /io/build-wheels.sh /io 6 7 8
