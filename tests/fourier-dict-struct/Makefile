INC = -I../..
SRC = ../../cm_dict_gen.c ./fourier-dict-struct.c
HDR = ../../cm_dict_gen.h
BIN = ./fourier-dict-struct
DEF = -DCM_DEBUG

$(BIN) : ${SRC} ${HDR}
	${CC} ${DEF} ${INC} ${SRC} -o ${BIN}

clean:
	rm ${BIN}
