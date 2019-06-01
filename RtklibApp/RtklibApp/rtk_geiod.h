#ifndef _RTKGEIOD_H_
#define _RTKGEIOD_H_

#include "rtk_common.h"

#ifdef __cplusplus
extern "C" {
#endif


/* geiod models --------------------------------------------------------------*/
extern int    opengeoid(int model, const char *file);
extern void   closegeoid(void);
extern double geoidh(const double *pos);



#ifdef __cplusplus
}
#endif
#endif 