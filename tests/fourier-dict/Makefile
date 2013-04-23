INC = -I../..
SRC = ../../cm_dict_gen.c ./fourier-dict-test.c
HDR = ../../cm_dict_gen.h
BIN = ./fourier-dict-test
DEF = -DCM_DEBUG

$(BIN) : ${SRC} ${HDR}
	${CC} ${DEF} ${INC} ${SRC} -o ${BIN}

clean:
	rm ${BIN}
