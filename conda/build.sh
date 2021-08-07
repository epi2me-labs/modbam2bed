#!/bin/bash

NAME=modbam2bed

export HTS_CONF_ARGS="--prefix=${PREFIX} --enable-libcurl --with-libdeflate --enable-plugins --enable-gcs --enable-s3"
export EXTRA_CFLAGS="-I$PREFIX/include"
export EXTRA_LDFLAGS="-L$PREFIX/lib"
export EXTRA_LIBS="-ldl -ldeflate"
make clean $NAME

mkdir -p $PREFIX/bin
cp $NAME $PREFIX/bin && chmod +x $PREFIX/bin/$NAME
