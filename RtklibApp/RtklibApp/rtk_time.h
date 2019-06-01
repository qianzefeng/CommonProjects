#ifndef _RTKTIME_H_
#define _RTKTIME_H_

#include "rtk_common.h"

#ifdef __cplusplus
extern "C" {
#endif


/* time and string functions -------------------------------------------------*/
extern double  str2num(const char *s, int i, int n);
extern int     str2time(const char *s, int i, int n, gtime_t *t);
extern void    time2str(gtime_t t, char *str, int n);
extern gtime_t epoch2time(const double *ep);
extern void    time2epoch(gtime_t t, double *ep);
extern gtime_t gpst2time(int week, double sec);
extern double  time2gpst(gtime_t t, int *week);
extern gtime_t gst2time(int week, double sec);
extern double  time2gst(gtime_t t, int *week);
extern gtime_t bdt2time(int week, double sec);
extern double  time2bdt(gtime_t t, int *week);
extern char    *time_str(gtime_t t, int n);

extern gtime_t timeadd  (gtime_t t, double sec);
extern double  timediff (gtime_t t1, gtime_t t2);
extern gtime_t gpst2utc (gtime_t t);
extern gtime_t utc2gpst (gtime_t t);
extern gtime_t gpst2bdt (gtime_t t);
extern gtime_t bdt2gpst (gtime_t t);
extern gtime_t timeget  (void);
extern void    timeset  (gtime_t t);
extern double  time2doy (gtime_t t);
extern double  utc2gmst (gtime_t t, double ut1_utc);
extern int     read_leaps(const char *file);

extern int adjgpsweek(int week);
extern unsigned int tickget(void);
extern void sleepms(int ms);

extern int reppath(const char *path, char *rpath, gtime_t time, const char *rov,
				   const char *base);
extern int reppaths(const char *path, char *rpaths[], int nmax, gtime_t ts,
					gtime_t te, const char *rov, const char *base);
	


#ifdef __cplusplus
}
#endif
#endif