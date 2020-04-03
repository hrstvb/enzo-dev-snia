/***********************************************************************
 /
 /  EXCEPTION CLASS
 /
 /  written by: Matthew Turk
 /  date:       June, 2009
 /  modified:
 /
 /  PURPOSE:
 /
 ************************************************************************/
#ifndef __ENZO_EXCEPTIONS__
#define __ENZO_EXCEPTIONS__

#include <execinfo.h>
#include <stdio.h>

#include <exception>
#include "ErrorExceptions.h"

// If we are using the new problem type initializers, we need to include these
// in a file that we know will be included before macros_and_parameters.h.
#ifdef NEW_PROBLEM_TYPES
#include <stdlib.h>
#include <stdio.h>
#include <map>
#include <string>
#endif

// Example from
//
// http://www.ibm.com/developerworks/linux/library/l-cppexcep.html
//

// This must be included BEFORE macros_and_parameters.h
// so we use int here

extern char current_error[255];

class EnzoFatalException
{
public:
	EnzoFatalException(const char *filename, const int line_number, const char *format, va_list varargs)
	{
		void * array[25];
		int nSize = backtrace(array, 25);
		char **symbols = backtrace_symbols(array, nSize);
		fprintf(stderr, "Caught fatal exception:\n\n   '");
		vfprintf(stderr, format, varargs);
		if(filename != NULL)
			fprintf(stderr, "'\n   at %s:%d\n\n", filename, line_number);
		else
			fprintf(stderr, "'\n");
		fprintf(stderr, "Backtrace:\n\n");

		for(int i = 0; i < nSize; i++)
		{
			fprintf(stderr, "BT symbol: %s\n", symbols[i]);
		}

		delete[] symbols;
	}

	EnzoFatalException(const char *error_msg, const char *filename = NULL, int line_number = 0)
	{
		void * array[25];
		int nSize = backtrace(array, 25);
		char ** symbols = backtrace_symbols(array, nSize);
		fprintf(stderr, "Caught fatal exception:\n\n");
		fprintf(stderr, "   '%s'\n", error_msg);
		if(filename != NULL)
			fprintf(stderr, "at %s:%d\n\n", filename, line_number);
		fprintf(stderr, "Backtrace:\n\n");

		for(int i = 0; i < nSize; i++)
		{
			fprintf(stderr, "BT symbol: %s\n", symbols[i]);
		}

		delete[] symbols;
	}

	void WriteDebuggingOutput()
	{

	}
};

/* We have two possibilities for throwing an exception.
 You can throw ENZO_FAIL, which takes no arguments and is safe to use with a
 semicolon, or you can throw ENZO_VFAIL which must NOT have a semicolon and
 comes enclosed in brackets.  You can supply format strings to ENZO_VFAIL.  */
#ifdef CONFIG_THROW_ABORT
#define ENZO_FAIL(A) raise(SIGABRT);
#define ENZO_VFAIL(A, ...) raise(SIGABRT);
#else
//#define ENZO_FAIL(A) throw(EnzoFatalException(A, __FILE__, __LINE__));
#define ENZO_THROW(A) do{throw(EnzoFatalException(__FILE__, __LINE__, A));}while(0)
////#define ENZO_VFAIL(format, ...) {snprintf(current_error, 254, format, ##__VA_ARGS__); throw(EnzoFatalException(current_error, __FILE__, __LINE__));}
#define ENZO_VTHROW(format, ...) do{snprintf(current_error, 254, format, ##__VA_ARGS__); throw(EnzoFatalException(current_error, __FILE__, __LINE__));}while(0)
#endif

#endif
