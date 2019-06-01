#ifndef _RTK_DATUM_H_
#define _RTK_DATUM_H_

#include "rtk_common.h"


#ifdef __cplusplus
extern "C" {
#endif


/* datum transformation ------------------------------------------------------*/
extern int loaddatump(const char *file);
extern int tokyo2jgd(double *pos);
extern int jgd2tokyo(double *pos);



#ifdef __cplusplus
}
#endif
#endif