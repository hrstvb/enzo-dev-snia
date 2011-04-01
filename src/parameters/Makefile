SRC = parameters.o parameters_Param.o parameters_Parameters.o
PARSE = lex.yy.o parse.tab.o
OBJ = $(SRC) $(PARSE)
INCLUDE = -I ../enzo

# compile everything
all: test_Parameters test_Parse
doc: 
	doxygen

# tests
test_Parameters: test_Parameters.o $(OBJ)
	g++ -o test_Parameters $(INCLUDE) test_Parameters.o $(OBJ)
test_Parse:test_Parse.o $(OBJ)
	g++ -o test_Parse $(INCLUDE) test_Parse.o $(OBJ)

# object files
parameters.o: parameters.cpp
	g++ -c parameters.cpp $(INCLUDE) 
parameters_Param.o: parameters_Param.cpp
	g++ -c parameters_Param.cpp $(INCLUDE) 
parameters_Parameters.o: parameters_Parameters.cpp
	g++ -c parameters_Parameters.cpp $(INCLUDE) 

clean:
	rm -rf *.o *~ test_Parameters test_Parse test.in test.out html

lex.yy.c: parse.l
	flex -t parse.l > lex.yy.c
parse.tab.c: parse.y
	bison -d -t -o parse.tab.c parse.y
