#ifndef SRC_ENZO_MYLIMITERS_H_
#define SRC_ENZO_MYLIMITERS_H_

//#include <float.h>

#define MYLIM_MAX (10)
#define MYLIM_FLOAT double
//#define MYLIM_FORMAT ".15f"
//#define MYLIM_EPSILON (DBL_EPSILON)
//#define MYLIM_TOLL (MYLIM_EPSILON*16)

#include "macros_and_parameters.h"

//extern MYLIM_FLOAT MyLimiterX12Y12[10][4];

void limiter1nInit();
MYLIM_FLOAT invLimiter(MYLIM_FLOAT dP01, MYLIM_FLOAT slope, int *err);

#endif /* SRC_ENZO_MYLIMITERS_H_ */
