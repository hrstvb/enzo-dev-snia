#ifndef SRC_ENZO_DEBUGMACROS_H_
#define SRC_ENZO_DEBUGMACROS_H_

#include "stdio.h"

//MACRO2STRING_ -- Not intended to be used directly.
#define MACRO2STRING_(x) #x

/*
 * MACRO2STR(X)
 * 		When X is another macro, e.g. #define X Y,
 * 		MACROTOSTR is replaced by "Y".  This is Y
 * 		enclosed in double quotes and any double
 * 		quotes in Y are eascaped.
 */
#define MACRO2STR(x) MACRO2STRING_(x)

// LINE2STR -- The current line number in double quotes.
#define LINE2STR MACRO2STR(__LINE__)

// FLN2STR -- Formatted string of the current file and line.
#define FLN2STR __FILE__ ":" LINE2STR

// MACRO_FMT -- Formated string of a macro name and replacement.
#define MACRO_FMT(macro) FLN2STR ": " #macro "=" MACRO2STR(macro)

/*
 * PRAGMA_MACRO(M) -- Can be used in a #pragma directive, e.g.
 * 		#pragma PRAGMA_MACRO(ANOTHER_MACRO)
 */
#define PRAGMA_MACRO(macro) "message(" MACRO_FMT(macro) ")"

/*
 * TRACE
 * TRACEF
 *		Enzo-specific macros providing an easy way to insert debug
 *		prints with the current file/line and MPI task number.
 *		TRACEF can print additional text or data.
 */

#ifndef ROOT_PROCESSOR
#define ROOT_PROCESSOR 0
#endif

extern long long MyProcessorNumber;

#define TRACE do{fprintf(stderr, "TRACE:(#%lld)" FLN2STR "\n", MyProcessorNumber);}while(0)
#define TRACEF(format, ...) do{fprintf(stderr, "TRACE:(#%lld)" FLN2STR ":" format "\n", MyProcessorNumber, ##__VA_ARGS__);}while(0)

#define TRACEG do{fprintf(stderr, "TRACE:(#%lld)g%lld(%p)::" FLN2STR "\n", MyProcessorNumber, this->ID, this);}while(0)
#define TRACEGF(format, ...) do{fprintf(stderr, "TRACE:(#%lld)g%lld(%p)::" FLN2STR ":" format "\n", MyProcessorNumber, this->ID, this, ##__VA_ARGS__);}while(0)

#define TRACE0 do{if(MyProcessorNumber==ROOT_PROCESSOR){TRACE;}}while(0)
#define TRACEF0(format, ...) do{if(MyProcessorNumber==ROOT_PROCESSOR){TRACEF(format, ##__VA_ARGS__);}}while(0)

#define TRACEG0 do{if(MyProcessorNumber==ROOT_PROCESSOR){TRACEG;}}while(0)
#define TRACEGF0(format, ...) do{if(MyProcessorNumber==ROOT_PROCESSOR){TRACEGF(format, ##__VA_ARGS__);}}while(0)

#define EXTRAFUNC(grid, i, ...) do{grid->ExtraFunction(FLN2STR, i, ##__VA_ARGS__);}while(0)
#define EXTRAFUNC0(grid, i, ...) do{if(MyProcessorNumber==ROOT_PROCESSOR){grid->ExtraFunction(FLN2STR, i, ##__VA_ARGS__);}}while(0)
#define EXTRAFUNCG(i, ...) do{ExtraFunction(FLN2STR, i, ##__VA_ARGS__);}while(0)
#define EXTRAFUNCG0(i, ...) do{if(MyProcessorNumber==ROOT_PROCESSOR){ExtraFunction(FLN2STR, i, ##__VA_ARGS__);}}while(0)

#endif /* SRC_ENZO_DEBUGMACROS_H_ */
