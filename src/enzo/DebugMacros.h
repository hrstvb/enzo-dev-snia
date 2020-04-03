#ifndef SRC_ENZO_DEBUGMACROS_H_
#define SRC_ENZO_DEBUGMACROS_H_

#define MACRO2STRING_(x) #x
#define MACRO2STR(x) MACRO2STRING_(x)
#define FILE2STR MACRO2STR(__FILE__)
#define LINE2STR MACRO2STR(__LINE__)
#define FLN2STR FILE2STR ":" LINE2STR
#define MACRO_FMT(macro) FLN2STR ": " #macro "=" MACRO2STR(macro)
#define PRAGMA_MACRO(macro) "message(" MACRO_FMT(macro) ")"

#define TRACE fprintf(stderr, "TRACE:(#%lld)" FLN2STR "\n", MyProcessorNumber);
#define TRACEF(format, ...) fprintf(stderr, "TRACE:(#%lld)" FLN2STR ":" format "\n", MyProcessorNumber, ##__VA_ARGS__);

#endif /* SRC_ENZO_DEBUGMACROS_H_ */
