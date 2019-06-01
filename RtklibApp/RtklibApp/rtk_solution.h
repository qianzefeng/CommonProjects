#ifndef _RTKSOLUTION_H_
#define _RTKSOLUTION_H_

#include "rtk_common.h"

#ifdef __cplusplus
extern "C" {
#endif

/* solution functions --------------------------------------------------------*/

extern void initsolbuf(solbuf_t *solbuf, int cyclic, int nmax);
extern void freesolbuf(solbuf_t *solbuf);
extern void freesolstatbuf(solstatbuf_t *solstatbuf);
extern sol_t *getsol(solbuf_t *solbuf, int index);
extern int addsol(solbuf_t *solbuf, const sol_t *sol);
extern int readsol (char *files[], int nfile, solbuf_t *sol);
extern int readsolt(char *files[], int nfile, gtime_t ts, gtime_t te, double tint, int qflag, solbuf_t *sol);
extern int readsolstat(char *files[], int nfile, solstatbuf_t *statbuf);
extern int readsolstatt(char *files[], int nfile, gtime_t ts, gtime_t te, double tint, solstatbuf_t *statbuf);
extern int inputsol(unsigned char data, gtime_t ts, gtime_t te, double tint, int qflag, const solopt_t *opt, solbuf_t *solbuf);
extern int outprcopts(unsigned char *buff, const prcopt_t *opt);
extern int outsolheads(unsigned char *buff, const solopt_t *opt);
extern int outsols  (unsigned char *buff, const sol_t *sol, const double *rb, const solopt_t *opt);
extern int outsolexs(unsigned char *buff, const sol_t *sol, const ssat_t *ssat, const solopt_t *opt);
extern void outprcopt(FILE *fp, const prcopt_t *opt);
extern void outsolhead(FILE *fp, const solopt_t *opt);
extern void outsol  (FILE *fp, const sol_t *sol, const double *rb, const solopt_t *opt);
extern void outsolex(FILE *fp, const sol_t *sol, const ssat_t *ssat, const solopt_t *opt);

extern int outnmea_rmc(unsigned char *buff, const sol_t *sol);
extern int outnmea_gga(unsigned char *buff, const sol_t *sol);
extern int outnmea_gsa(unsigned char *buff, const sol_t *sol, const ssat_t *ssat);
extern int outnmea_gsv(unsigned char *buff, const sol_t *sol, const ssat_t *ssat);

#ifdef __cplusplus
}
#endif
#endif 