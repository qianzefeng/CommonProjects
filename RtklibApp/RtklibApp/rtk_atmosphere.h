#ifndef _RTKATMOSPHERE_H_
#define _RTKATMOSPHERE_H_

#include "rtk_common.h"

#ifdef __cplusplus
extern "C" {
#endif

/* atmosphere models ---------------------------------------------------------*/
extern double ionmodel(gtime_t t, const double *ion, const double *pos, const double *azel);
extern double ionmapf(const double *pos, const double *azel);
extern double ionppp(const double *pos, const double *azel, double re, double hion, double *pppos);
extern double tropmodel(gtime_t time, const double *pos, const double *azel, double humi);
extern double tropmapf(gtime_t time, const double *pos, const double *azel, double *mapfw);
extern int iontec(gtime_t time, const nav_t *nav, const double *pos, const double *azel, int opt, double *delay, double *var);
extern void readtec(const char *file, nav_t *nav, int opt);
extern int ionocorr(gtime_t time, const nav_t *nav, int sat, const double *pos, const double *azel, int ionoopt, double *ion, double *var);
extern int tropcorr(gtime_t time, const nav_t *nav, const double *pos, const double *azel, int tropopt, double *trp, double *var);

#ifdef __cplusplus
}
#endif
#endif 