#ifndef _RTKKML_H_
#define _RTKKML_H_

#include "rtk_common.h"

#ifdef __cplusplus
extern "C" {
#endif

/* google earth kml converter ------------------------------------------------*/
extern int convkml(const char *infile, const char *outfile, gtime_t ts,
				   gtime_t te, double tint, int qflg, double *offset,
				   int tcolor, int pcolor, int outalt, int outtime);

#ifdef __cplusplus
}
#endif
#endif 