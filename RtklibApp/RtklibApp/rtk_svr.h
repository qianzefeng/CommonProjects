#ifndef _RTKSVR_H_
#define _RTKSVR_H_

#include "rtk_common.h"

#ifdef __cplusplus
extern "C" {
#endif


/* rtk server functions ------------------------------------------------------*/
extern int  rtksvrinit  (rtksvr_t *svr);
extern void rtksvrfree  (rtksvr_t *svr);
extern int  rtksvrstart (rtksvr_t *svr, int cycle, int buffsize, int *strs,
                         char **paths, int *formats, int navsel, char **cmds,
                         char **rcvopts, int nmeacycle, int nmeareq,
                         const double *nmeapos, prcopt_t *prcopt,
                         solopt_t *solopt, stream_t *moni);
extern void rtksvrstop  (rtksvr_t *svr, char **cmds);
extern int  rtksvropenstr(rtksvr_t *svr, int index, int str, const char *path,
                          const solopt_t *solopt);
extern void rtksvrclosestr(rtksvr_t *svr, int index);
extern void rtksvrlock  (rtksvr_t *svr);
extern void rtksvrunlock(rtksvr_t *svr);
extern int  rtksvrostat (rtksvr_t *svr, int type, gtime_t *time, int *sat,
                         double *az, double *el, int **snr, int *vsat);
extern void rtksvrsstat (rtksvr_t *svr, int *sstat, char *msg);



#ifdef __cplusplus
}
#endif
#endif 