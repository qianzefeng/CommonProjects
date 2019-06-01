#ifndef _RTK_MAT_H_
#define _RTK_MAT_H_

#include "rtk_common.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef void fatalfunc_t(const char *);		/* fatal callback function type */
extern void add_fatal(fatalfunc_t *func);

/* matrix and vector functions -----------------------------------------------*/
extern double *mat  (int n, int m);
extern int    *imat (int n, int m);
extern double *zeros(int n, int m);
extern double *eye  (int n);
extern double dot (const double *a, const double *b, int n);
extern double norm(const double *a, int n);
extern void cross3(const double *a, const double *b, double *c);
extern int  normv3(const double *a, double *b);
extern void matcpy(double *A, const double *B, int n, int m);
extern void matmul(const char *tr, int n, int k, int m, double alpha,
				   const double *A, const double *B, double beta, double *C);
extern int  matinv(double *A, int n);
extern int  solve (const char *tr, const double *A, const double *Y, int n,
				   int m, double *X);
extern int  lsq   (const double *A, const double *y, int n, int m, double *x,
				   double *Q);
extern int  filter(double *x, double *P, const double *H, const double *v,
				   const double *R, int n, int m);
extern int  smoother(const double *xf, const double *Qf, const double *xb,
					 const double *Qb, int n, double *xs, double *Qs);
extern void matprint (const double *A, int n, int m, int p, int q);
extern void matfprint(const double *A, int n, int m, int p, int q, FILE *fp);


#ifdef __cplusplus
}
#endif
#endif