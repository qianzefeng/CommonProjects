#ifndef _RTKTRACE_H_
#define _RTKTRACE_H_

#include "rtk_common.h"

//#define TRACE

#ifdef __cplusplus
extern "C" {
#endif

/* debug trace functions -----------------------------------------------------*/
extern void traceopen(const char *file);
extern void traceclose(void);
extern void tracelevel(int level);
extern void trace    (int level, const char *format, ...);
extern void tracet   (int level, const char *format, ...);
extern void tracemat (int level, const double *A, int n, int m, int p, int q);
extern void traceobs (int level, const obsd_t *obs, int n);
extern void tracenav (int level, const nav_t *nav);
extern void tracegnav(int level, const nav_t *nav);
extern void tracehnav(int level, const nav_t *nav);
extern void tracepeph(int level, const nav_t *nav);
extern void tracepclk(int level, const nav_t *nav);
extern void traceb   (int level, const unsigned char *p, int n);


/* application defined functions ---------------------------------------------*/
extern int showmsg(char *format,...);
extern void settspan(gtime_t ts, gtime_t te);
extern void settime(gtime_t time);

#ifdef __cplusplus
}
#endif
#endif 