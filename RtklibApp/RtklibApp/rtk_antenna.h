#ifndef _RTKANTENNA_H_
#define _RTKANTENNA_H_

#include "rtk_common.h"

#ifdef __cplusplus
extern "C" {
#endif


/* antenna models ------------------------------------------------------------*/
extern int  readpcv(const char *file, pcvs_t *pcvs);
extern pcv_t *searchpcv(int sat, const char *type, gtime_t time, const pcvs_t *pcvs);
extern void antmodel(const pcv_t *pcv, const double *del, const double *azel, int opt, double *dant);
extern void antmodel_s(const pcv_t *pcv, double nadir, double *dant);


#ifdef __cplusplus
}
#endif
#endif 