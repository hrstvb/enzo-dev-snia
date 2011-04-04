SRC_CXX  = parameters.cpp parameters_Param.cpp parameters_Parameters.cpp
SRC_C    = parse.tab.c lex.yy.c

OBJ  = $(patsubst %.cpp,%.o,$(SRC_CXX)) $(patsubst %.c,%.o,$(SRC_C))
EXE    = test_Parse test_Parameters


all: run-exe

#=======================================================================
# TARGETS
#=======================================================================

.PHONY: compile
compile: $(EXE)

.PHONY: run-exe
run-exe: $(EXE)
	-./test_Parse < input/parse_groups.in > output/parse_groups.out
	-./test_Parse < input/parse_implosion.in > output/parse_implosion.out
	-./test_Parse < input/parse_integer.in > output/parse_integer.out
	-./test_Parse < input/parse_list.in > output/parse_list.out
	-./test_Parse < input/parse_logical.in > output/parse_logical.out
	-./test_Parse < input/parse_scalar.in > output/parse_scalar.out
	test_Parameters > output/test_Parameters.out
	-egrep "FAIL|pass" output/*

.PHONY: doc
doc:
	doxygen  # creates html/index.html

#-----------------------------------------------------------------------

test_Parameters: test_Parameters.o $(OBJ)
	$(CXX) -o test_Parameters test_Parameters.o $(OBJ)

test_Parse: test_Parse.o $(OBJ)
	$(CXX) -o test_Parse test_Parse.o $(OBJ)

parse.tab.o:parse.tab.c
	$(CC) parse.tab.c -c
lex.yy.o:lex.yy.c
	$(CC) lex.yy.c -c

lex.yy.c:parse.l
	flex parse.l
parse.tab.h: parse.y
	bison -d parse.y 
parse.tab.c: parse.y
	bison -d parse.y 

.PHONY: clean
clean:
	-@(rm -rf output/* *~ *.o test.in test.out html $(EXE))


