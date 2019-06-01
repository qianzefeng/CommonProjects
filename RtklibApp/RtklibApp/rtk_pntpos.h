#ifndef _RTK_PNTPOS_H_
#define _RTK_PNTPOS_H_

#include "rtk_common.h"

#ifdef __cplusplus
extern "C" {
#endif

/* standard positioning ------------------------------------------------------*/
extern int pntpos(const obsd_t *obs, int n, const nav_t *nav,
				  const prcopt_t *opt, sol_t *sol, double *azel,
				  ssat_t *ssat, char *msg);

#ifdef __cplusplus
}
#endif
#endif