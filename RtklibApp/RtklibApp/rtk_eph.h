#ifndef _RTKEPH_H_
#define _RTKEPH_H_

#include "rtk_common.h"

#ifdef __cplusplus
extern "C" {
#endif


/* ephemeris and clock functions ---------------------------------------------*/
/* precise ephemeris and clock functions--------------------------------------*/


extern double eph2clk (gtime_t time, const eph_t  *eph);
extern double geph2clk(gtime_t time, const geph_t *geph);
extern double seph2clk(gtime_t time, const seph_t *seph);
extern void eph2pos (gtime_t time, const eph_t  *eph,  double *rs, double *dts,
                     double *var);
extern void geph2pos(gtime_t time, const geph_t *geph, double *rs, double *dts,
                     double *var);
extern void seph2pos(gtime_t time, const seph_t *seph, double *rs, double *dts,
                     double *var);
extern int  peph2pos(gtime_t time, int sat, const nav_t *nav, int opt,
                     double *rs, double *dts, double *var);
extern void satantoff(gtime_t time, const double *rs, int sat, const nav_t *nav,
                      double *dant);
extern int  satpos(gtime_t time, gtime_t teph, int sat, int ephopt,
                   const nav_t *nav, double *rs, double *dts, double *var,
                   int *svh);
extern void satposs(gtime_t time, const obsd_t *obs, int n, const nav_t *nav,
                    int sateph, double *rs, double *dts, double *var, int *svh);
extern void readsp3(const char *file, nav_t *nav, int opt);
extern int  readsap(const char *file, gtime_t time, nav_t *nav);
extern int  readdcb(const char *file, nav_t *nav, const sta_t *sta);
extern int  readfcb(const char *file, nav_t *nav);
extern void alm2pos(gtime_t time, const alm_t *alm, double *rs, double *dts);

extern int tle_read(const char *file, tle_t *tle);
extern int tle_name_read(const char *file, tle_t *tle);
extern int tle_pos(gtime_t time, const char *name, const char *satno,
                   const char *desig, const tle_t *tle, const erp_t *erp,
                   double *rs);


#ifdef __cplusplus
}
#endif
#endif 