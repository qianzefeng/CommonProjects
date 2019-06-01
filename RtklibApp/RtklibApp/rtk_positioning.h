#ifndef _RTKPOSITIONING_H_
#define _RTKPOSITIONING_H_

#include "rtk_common.h"

#ifdef __cplusplus
extern "C" {
#endif


/* positioning models --------------------------------------------------------*/
extern double satwavelen(int sat, int frq, const nav_t *nav);
extern double satazel(const double *pos, const double *e, double *azel);
extern double geodist(const double *rs, const double *rr, double *e);
extern void dops(int ns, const double *azel, double elmin, double *dop);
extern void csmooth(obs_t *obs, int ns);

#ifdef __cplusplus
}
#endif
#endif 