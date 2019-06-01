#ifndef _RTKLAMBDA_H_
#define _RTKLAMBDA_H_

#include "rtk_common.h"

#ifdef __cplusplus
extern "C" {
#endif

/* integer ambiguity resolution ----------------------------------------------*/
extern int lambda(int n, int m, const double *a, const double *Q, double *F, double *s);
extern int lambda_reduction(int n, const double *Q, double *Z);
extern int lambda_search(int n, int m, const double *a, const double *Q, double *F, double *s);


#ifdef __cplusplus
}
#endif
#endif