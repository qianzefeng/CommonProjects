#ifndef _RTKSTREAMSVR_H_
#define _RTKSTREAMSVR_H_

#include "rtk_common.h"

#ifdef __cplusplus
extern "C" {
#endif


/* stream server functions ---------------------------------------------------*/
extern void strsvrinit (strsvr_t *svr, int nout);
extern int  strsvrstart(strsvr_t *svr, int *opts, int *strs, char **paths,
						strconv_t **conv, const char *cmd,
						const double *nmeapos);
extern void strsvrstop (strsvr_t *svr, const char *cmd);
extern void strsvrstat (strsvr_t *svr, int *stat, int *byte, int *bps, char *msg);
extern strconv_t *strconvnew(int itype, int otype, const char *msgs, int staid,
							 int stasel, const char *opt);
extern void strconvfree(strconv_t *conv);


#ifdef __cplusplus
}
#endif
#endif 