#ifndef _RTKRTCM_H_
#define _RTKRTCM_H_

#include "rtk_common.h"

#ifdef __cplusplus
extern "C" {
#endif

extern int decode_rtcm2(rtcm_t *rtcm);
extern int decode_rtcm3(rtcm_t *rtcm);
extern int encode_rtcm3(rtcm_t *rtcm, int type, int sync);

/* rtcm functions ------------------------------------------------------------*/
extern int init_rtcm   (rtcm_t *rtcm);
extern void free_rtcm  (rtcm_t *rtcm);
extern int input_rtcm2 (rtcm_t *rtcm, unsigned char data);
extern int input_rtcm3 (rtcm_t *rtcm, unsigned char data);
extern int input_rtcm2f(rtcm_t *rtcm, FILE *fp);
extern int input_rtcm3f(rtcm_t *rtcm, FILE *fp);
extern int gen_rtcm2   (rtcm_t *rtcm, int type, int sync);
extern int gen_rtcm3   (rtcm_t *rtcm, int type, int sync);



#ifdef __cplusplus
}
#endif
#endif 