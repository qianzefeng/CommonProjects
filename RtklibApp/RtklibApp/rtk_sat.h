#ifndef _RTKSAT_H_
#define _RTKSAT_H_

#include "rtk_common.h"

#ifdef __cplusplus
extern "C" {
#endif


/* satellites, systems, codes functions --------------------------------------*/
extern int  satno   (int sys, int prn);
extern int  satsys  (int sat, int *prn);
extern int  satid2no(const char *id);
extern void satno2id(int sat, char *id);
extern unsigned char obs2code(const char *obs, int *freq);
extern char *code2obs(unsigned char code, int *freq);
extern int  satexclude(int sat, int svh, const prcopt_t *opt);
extern int  testsnr(int base, int freq, double el, double snr, const snrmask_t *mask);
extern void setcodepri(int sys, int freq, const char *pri);
extern int  getcodepri(int sys, unsigned char code, const char *opt);


#ifdef __cplusplus
}
#endif
#endif 