#!/bin/bash
gcc -I/usr/local/include -I../.. /usr/local/lib/libsndfile.dylib  \
    /usr/local/lib/libglpk.dylib \
    /usr/local/lib/libgsl.dylib \
    /usr/local/lib/libgslcblas.dylib \
    -I/Users/nicklocal/Documents/mumt-622/code/libcm/inc/ \
    /Users/nicklocal/Documents/mumt-622/code/libcm/sndfio.c \
    -D_CM_DEBUG \
    ../../cm_dict_gen.c ./cm-mp-bp-cmp.c ../../cm_bp.c \
    ../../cm_mp.c -o ./cm-mp-bp-cmp
