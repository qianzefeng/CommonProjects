#ifndef _RTKQZSS_H_
#define _RTKQZSS_H_

#include "rtk_common.h"

#ifdef __cplusplus
extern "C" {
#endif

/* qzss lex functions --------------------------------------------------------*/
extern int lexupdatecorr(const lexmsg_t *msg, nav_t *nav, gtime_t *tof);
extern int lexreadmsg(const char *file, int sel, lex_t *lex);
extern void lexoutmsg(FILE *fp, const lexmsg_t *msg);
extern int lexconvbin(int type, int format, const char *infile,
					  const char *outfile);
extern int lexeph2pos(gtime_t time, int sat, const nav_t *nav, double *rs,
					  double *dts, double *var);
extern int lexioncorr(gtime_t time, const nav_t *nav, const double *pos,
					  const double *azel, double *delay, double *var);


#ifdef __cplusplus
}
#endif
#endif 