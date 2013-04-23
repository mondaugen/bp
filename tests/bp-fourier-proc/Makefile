INC = -I../.. -I/Users/nicklocal/Documents/mumt-622/code/libcm/inc/\
	  -L/usr/local/lib/ -I/usr/local/include/ -I/usr/include
SRC = ../../cm_dict_gen.c ./bp-fourier-proc.c\
	  /Users/nicholasesterer/Documents/code/c/computermusic/libcm/src/sndfio.c
HDR = ../../cm_dict_gen.h\
	  /Users/nicholasesterer/Documents/code/c/computermusic/libcm/src/inc/cm/sndfio.h
LIB = /usr/lib/libsndfile.dylib /usr/local/lib/libglpk.dylib
BIN = ./bp-fourier-proc
DEF = #-DCM_DEBUG

$(BIN) : ${SRC} ${HDR}
	${CC} ${DEF} ${INC} ${SRC} ${LIB} -o ${BIN}

clean:
	rm ${BIN}
