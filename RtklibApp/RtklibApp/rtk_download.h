#ifndef _RTKANTENNA_H_
#define _RTKANTENNA_H_

#include "rtk_common.h"

#ifdef __cplusplus
extern "C" {
#endif

/* downloader functions ------------------------------------------------------*/
extern int dl_readurls(const char *file, char **types, int ntype, url_t *urls,
					   int nmax);
extern int dl_readstas(const char *file, char **stas, int nmax);
extern int dl_exec(gtime_t ts, gtime_t te, double ti, int seqnos, int seqnoe,
				   const url_t *urls, int nurl, char **stas, int nsta,
				   const char *dir, const char *usr, const char *pwd,
				   const char *proxy, int opts, char *msg, FILE *fp);
extern void dl_test(gtime_t ts, gtime_t te, double ti, const url_t *urls,
					int nurl, char **stas, int nsta, const char *dir,
					int ncol, int datefmt, FILE *fp);


#ifdef __cplusplus
}
#endif
#endif 