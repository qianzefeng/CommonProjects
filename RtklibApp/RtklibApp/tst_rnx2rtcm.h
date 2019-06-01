#ifndef _RTKRNX2RTCM_H_
#define _RTKRNX2RTCM_H_

#include "rtk_common.h"

#ifdef __cplusplus
extern "C" {
#endif


extern int conv_rtcm(const int *type, int n, const char *outfile,
					 const obs_t *obs, const nav_t *nav, const sta_t *sta,
					 int staid);



#ifdef __cplusplus
}
#endif
#endif