#ifndef _RTKSBAS_H_
#define _RTKSBAS_H_

#include "rtk_common.h"

#ifdef __cplusplus
extern "C" {
#endif


/* sbas functions ------------------------------------------------------------*/
extern int  sbsreadmsg (const char *file, int sel, sbs_t *sbs);
extern int  sbsreadmsgt(const char *file, int sel, gtime_t ts, gtime_t te,
						sbs_t *sbs);
extern void sbsoutmsg(FILE *fp, sbsmsg_t *sbsmsg);
extern int  sbsdecodemsg(gtime_t time, int prn, const unsigned int *words,
						 sbsmsg_t *sbsmsg);
extern int sbsupdatecorr(const sbsmsg_t *msg, nav_t *nav);
extern int sbssatcorr(gtime_t time, int sat, const nav_t *nav, double *rs,
					  double *dts, double *var);
extern int sbsioncorr(gtime_t time, const nav_t *nav, const double *pos,
					  const double *azel, double *delay, double *var);
extern double sbstropcorr(gtime_t time, const double *pos, const double *azel,
						  double *var);




#ifdef __cplusplus
}
#endif
#endif 