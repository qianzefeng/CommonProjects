#ifndef _RTK_POSTPOS_H_
#define _RTK_POSTPOS_H_

#include "rtk_common.h"

#ifdef __cplusplus
extern "C" {
#endif

/* post-processing positioning -----------------------------------------------*/
extern int postpos(gtime_t ts, gtime_t te, double ti, double tu,
				   const prcopt_t *popt, const solopt_t *sopt,
				   const filopt_t *fopt, char **infile, int n, char *outfile,
				   const char *rov, const char *base);


#ifdef __cplusplus
}
#endif
#endif