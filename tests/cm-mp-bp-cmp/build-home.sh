#!/bin/bash
gcc -I/usr/local/include -I../.. /usr/lib/libsndfile.dylib  \
    /usr/local/lib/libglpk.dylib \
    /usr/local/lib/libgsl.dylib \
    /usr/local/lib/libgslcblas.dylib \
    -I/Users/nicholasesterer/Documents/code/c/computermusic/libcm/src/inc/ \
    -D_CM_DEBUG \
    /Users/nicholasesterer/Documents/code/c/computermusic/libcm/src/sndfio.c \
    ../../cm_dict_gen.c ./cm-mp-bp-cmp.c ../../cm_bp.c \
    ../../cm_mp.c -o ./cm-mp-bp-cmp
