#!/bin/bash
gcc -I/usr/local/include -I../.. /usr/local/lib/libsndfile.dylib  \
    /usr/local/lib/libglpk.dylib \
    -I/Users/nicklocal/Documents/mumt-622/code/libcm/inc/ \
    /Users/nicklocal/Documents/mumt-622/code/libcm/sndfio.c \
    -D_CM_DEBUG \
    ../../cm_dict_gen.c ./cm-bp-test.c ../../cm_bp.c -o ./cm-bp-test
