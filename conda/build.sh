#!/bin/bash

NAME=modbam2bed

## self-built htslib
#export HTS_CONF_ARGS="--prefix=${PREFIX} --enable-libcurl --with-libdeflate --enable-plugins --enable-gcs --enable-s3"
#export EXTRA_CFLAGS="-I$PREFIX/include"
#export EXTRA_LDFLAGS="-L$PREFIX/lib"
#export EXTRA_LIBS="-ldl -lhts -ldeflate"
##export STATIC_HTSLIB=""

# just link to htslib from bioconda
export EXTRA_CFLAGS="-I$PREFIX/include"
export STATIC_HTSLIB=""
export EXTRA_LDFLAGS="-L$PREFIX/lib"
export EXTRA_LIBS="-ldl -lhts"

OS=$(uname)
if [[ "$OS" == "Darwin" ]]; then
    echo "Setting Darwin args"
    export ARGP=${PREFIX}/lib/libargp.a
    export EXTRA_CFLAGS="${EXTRA_CFLAGS} -isysroot ${CONDA_BUILD_SYSROOT} -mmacosx-version-min=${MACOSX_DEPLOYMENT_TARGET}"
fi

make clean $NAME

mkdir -p $PREFIX/bin
cp $NAME $PREFIX/bin && chmod +x $PREFIX/bin/$NAME
