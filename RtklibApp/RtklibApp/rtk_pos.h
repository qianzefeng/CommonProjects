#ifndef _RTKPOS_H_
#define _RTKPOS_H_

#include "rtk_common.h"

#ifdef __cplusplus
extern "C" {
#endif

/* precise positioning -------------------------------------------------------*/
extern void rtkinit(rtk_t *rtk, const prcopt_t *opt);
extern void rtkfree(rtk_t *rtk);
extern int  rtkpos (rtk_t *rtk, const obsd_t *obs, int nobs, const nav_t *nav);
extern int  rtkopenstat(const char *file, int level);
extern void rtkclosestat(void);
extern int  rtkoutstat(rtk_t *rtk, char *buff);


#ifdef __cplusplus
}
#endif
#endif 