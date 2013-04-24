# Makefile for pdlblp-test.c

INC = -I/usr/local/include -I../../ -I. -L/usr/local/lib
SRC = pdlblp-test.c ../../cm_pdlblp.c
OBJ = pdlblp-test.o cm_pdlblp.o
HDR = ../../cm_pdlblp.h
LIB = /usr/local/lib/libgsl.dylib /usr/local/lib/libgslcblas.dylib
#LIB = -lgsl -lgslcblas
BIN = pdlblp
OPT = -ggdb3
CC  = gcc

$(BIN) : $(OBJ)
	$(CC) $(OPT) $(INC) $(LIB) $(OBJ) -o $(BIN)

pdlblp-test.o : pdlblp-test.c
	$(CC) $(INC) $(OPT) -c $< -o $@

cm_pdlblp.o : ../../cm_pdlblp.c
	$(CC) $(INC) $(OPT) -c $< -o $@

clean :
	rm $(BIN)
