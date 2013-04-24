#!/bin/bash
gcc -I/usr/local/include -I../.. /usr/local/lib/libsndfile.dylib  \
    -I/Users/nicklocal/Documents/mumt-622/code/libcm/inc/ \
    /Users/nicklocal/Documents/mumt-622/code/libcm/sndfio.c \
    -D_CM_DEBUG \
    ../../cm_dict_gen.c ./dict-sndfile-loading.c -o ./dict-sndfile-loading
