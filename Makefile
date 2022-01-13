OS := $(shell uname)
ifeq ($(OS), Darwin)
	# mainly for dev builds using homebrew things
    EXTRA_LDFLAGS ?= -L/usr/local/Cellar/openssl@1.1/1.1.1k/lib
    ARGP ?= /usr/local/Cellar/argp-standalone/1.3/lib/libargp.a
endif

CC ?= gcc
CFLAGS ?= -fpic -msse3 -O3 -std=c99
STATIC_HTSLIB ?= htslib/libhts.a
EXTRA_CFLAGS ?=
EXTRA_LDFLAGS ?=
EXTRA_LIBS ?=
HTS_CONF_ARGS ?=
NOTHREADS ?=
ifeq ($(NOTHREADS), 1)
	CFLAGS += -DNOTHREADS
endif

VALGRIND ?= valgrind


.PHONY: default
default: modbam2bed


htslib/libhts.a:
	@echo Compiling $(@F)
	cd htslib/ \
		&& autoheader \
		&& autoconf \
		&& CFLAGS="$(CFLAGS) $(EXTRA_CFLAGS)" ./configure $(HTS_CONF_ARGS) \
		&& make -j 4


.PHONY: clean_htslib
clean_htslib:
	rm -rf htslib/autom4te.cache/ 
	cd htslib && make clean || exit 0


%.o: src/%.c
	mkdir -p obj && \
		$(CC) -c -pthread -Wall -fstack-protector-strong -D_FORTIFY_SOURCE=2 $(CFLAGS) \
		-Isrc -Ihtslib $(EXTRA_CFLAGS) $^ -o $@

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
		./modbam2bed -b 0.66 -a 0.33 -t 2 -r ecoli1 test_data/ecoli.fasta.gz test_data/400ecoli.bam > /dev/null


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
	${IN_VENV} && python setup.py develop

.PHONY: clean_python
clean_python: clean_obj
	rm -rf dist build modbampy.egg-info pymod.a libmodbampy.abi3.so ${VENV}

pymod.a: common.o bamiter.o counts.o args.o 
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
wheels: clean 
	docker run -v `pwd`:/io quay.io/pypa/manylinux2010_x86_64 /io/build-wheels.sh /io 6 7 8
