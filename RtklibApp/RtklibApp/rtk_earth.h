#ifndef _RTKEARTH_H_
#define _RTKEARTH_H_

#include "rtk_common.h"

#ifdef __cplusplus
extern "C" {
#endif


/* earth tide models ---------------------------------------------------------*/
extern void sunmoonpos(gtime_t tutc, const double *erpv, double *rsun,
					   double *rmoon, double *gmst);
extern void tidedisp(gtime_t tutc, const double *rr, int opt, const erp_t *erp,
					 const double *odisp, double *dr);

#ifdef __cplusplus
}
#endif
#endif 