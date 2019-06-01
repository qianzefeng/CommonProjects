#ifndef _RTKPPP_H_
#define _RTKPPP_H_

#include "rtk_common.h"

#ifdef __cplusplus
extern "C" {
#endif

/* precise point positioning -------------------------------------------------*/
extern void pppos(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav);
extern int pppnx(const prcopt_t *opt);
extern int pppoutstat(rtk_t *rtk, char *buff);

extern int ppp_ar(rtk_t *rtk, const obsd_t *obs, int n, int *exc,
				  const nav_t *nav, const double *azel, double *x, double *P);

extern int pppcorr_read(pppcorr_t *corr, const char *file);
extern void pppcorr_free(pppcorr_t *corr);
extern int pppcorr_trop(const pppcorr_t *corr, gtime_t time, const double *pos,
						double *ztd, double *std);
extern int pppcorr_stec(const pppcorr_t *corr, gtime_t time, const double *pos,
						double *ion, double *std);


#ifdef __cplusplus
}
#endif
#endif 