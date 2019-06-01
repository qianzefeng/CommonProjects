#ifndef _RTKPLATFORM_H_
#define _RTKPLATFORM_H_

#include "rtk_common.h"

#ifdef __cplusplus
extern "C" {
#endif

/* platform dependent functions ----------------------------------------------*/
extern int execcmd(const char *cmd);
extern int expath (const char *path, char *paths[], int nmax);
extern void createdir(const char *path);


#ifdef __cplusplus
}
#endif
#endif 