#ifndef _RTK_COORD_H_
#define _RTK_COORD_H_

#include "rtk_common.h"


#ifdef __cplusplus
extern "C" {
#endif

/* coordinates transformation ------------------------------------------------*/
extern void   deg2dms (double deg, double *dms);
extern double dms2deg(const double *dms);

extern void ecef2pos(const double *r, double *pos);
extern void pos2ecef(const double *pos, double *r);
extern void ecef2enu(const double *pos, const double *r, double *e);
extern void enu2ecef(const double *pos, const double *e, double *r);
extern void covenu  (const double *pos, const double *P, double *Q);
extern void covecef (const double *pos, const double *Q, double *P);
extern void xyz2enu (const double *pos, double *E);
extern void eci2ecef(gtime_t tutc, const double *erpv, double *U, double *gmst);

extern void ast_args(double t, double *f);

#ifdef __cplusplus
}
#endif
#endif