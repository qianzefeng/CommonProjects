#ifndef _RTKRINEX_H_
#define _RTKRINEX_H_

#include "rtk_common.h"


#ifdef __cplusplus
extern "C" {
#endif


/* rinex functions -----------------------------------------------------------*/
extern int readrnx (const char *file, int rcv, const char *opt, obs_t *obs, nav_t *nav, sta_t *sta);
extern int readrnxt(const char *file, int rcv, gtime_t ts, gtime_t te, double tint, const char *opt, obs_t *obs, nav_t *nav, sta_t *sta);
extern int readrnxc(const char *file, nav_t *nav);
extern int outrnxobsh(FILE *fp, const rnxopt_t *opt, const nav_t *nav);
extern int outrnxobsb(FILE *fp, const rnxopt_t *opt, const obsd_t *obs, int n, int epflag);
extern int outrnxnavh (FILE *fp, const rnxopt_t *opt, const nav_t *nav);
extern int outrnxgnavh(FILE *fp, const rnxopt_t *opt, const nav_t *nav);
extern int outrnxhnavh(FILE *fp, const rnxopt_t *opt, const nav_t *nav);
extern int outrnxlnavh(FILE *fp, const rnxopt_t *opt, const nav_t *nav);
extern int outrnxqnavh(FILE *fp, const rnxopt_t *opt, const nav_t *nav);
extern int outrnxcnavh(FILE *fp, const rnxopt_t *opt, const nav_t *nav);
extern int outrnxnavb (FILE *fp, const rnxopt_t *opt, const eph_t *eph);
extern int outrnxgnavb(FILE *fp, const rnxopt_t *opt, const geph_t *geph);
extern int outrnxhnavb(FILE *fp, const rnxopt_t *opt, const seph_t *seph);
extern int convrnx(int format, rnxopt_t *opt, const char *file, char **ofile);
extern int  init_rnxctr (rnxctr_t *rnx);
extern void free_rnxctr (rnxctr_t *rnx);
extern int  open_rnxctr (rnxctr_t *rnx, FILE *fp);
extern int  input_rnxctr(rnxctr_t *rnx, FILE *fp);

extern int uncompress(const char *file, char *uncfile);

#ifdef __cplusplus
}
#endif
#endif