#ifndef _RTK_IO_H_
#define _RTK_IO_H_

#include "rtk_common.h"

#ifdef __cplusplus
extern "C" {
#endif

/* input and output functions ------------------------------------------------*/
extern void readpos(const char *file, const char *rcv, double *pos);
extern int  sortobs(obs_t *obs);
extern void uniqnav(nav_t *nav);
extern int  screent(gtime_t time, gtime_t ts, gtime_t te, double tint);
extern int  readnav(const char *file, nav_t *nav);
extern int  savenav(const char *file, const nav_t *nav);
extern void freeobs(obs_t *obs);
extern void freenav(nav_t *nav, int opt);
extern int  readblq(const char *file, const char *sta, double *odisp);
extern int  readerp(const char *file, erp_t *erp);
extern int  geterp (const erp_t *erp, gtime_t time, double *val);


#ifdef __cplusplus
}
#endif
#endif