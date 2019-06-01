#ifndef _RTKRAW_H_
#define _RTKRAW_H_

#include "rtk_common.h"

#ifdef __cplusplus
extern "C" {
#endif


/* receiver raw data functions -----------------------------------------------*/
extern unsigned int getbitu(const unsigned char *buff, int pos, int len);
extern int          getbits(const unsigned char *buff, int pos, int len);
extern void setbitu(unsigned char *buff, int pos, int len, unsigned int data);
extern void setbits(unsigned char *buff, int pos, int len, int data);
extern unsigned int crc32  (const unsigned char *buff, int len);
extern unsigned int crc24q (const unsigned char *buff, int len);
extern unsigned short crc16(const unsigned char *buff, int len);

extern int decode_word (unsigned int word, unsigned char *data);
extern int decode_frame(const unsigned char *buff, eph_t *eph, alm_t *alm,
						double *ion, double *utc, int *leaps);
extern int test_glostr(const unsigned char *buff);
extern int decode_glostr(const unsigned char *buff, geph_t *geph);
extern int decode_bds_d1(const unsigned char *buff, eph_t *eph);
extern int decode_bds_d2(const unsigned char *buff, eph_t *eph);



extern int init_raw   (raw_t *raw);
extern void free_raw  (raw_t *raw);
extern int input_raw  (raw_t *raw, int format, unsigned char data);
extern int input_rawf (raw_t *raw, int format, FILE *fp);

extern int input_oem4  (raw_t *raw, unsigned char data);
extern int input_oem3  (raw_t *raw, unsigned char data);
extern int input_ubx   (raw_t *raw, unsigned char data);
extern int input_ss2   (raw_t *raw, unsigned char data);
extern int input_cres  (raw_t *raw, unsigned char data);
extern int input_stq   (raw_t *raw, unsigned char data);
extern int input_gw10  (raw_t *raw, unsigned char data);
extern int input_javad (raw_t *raw, unsigned char data);
extern int input_nvs   (raw_t *raw, unsigned char data);
extern int input_bnx   (raw_t *raw, unsigned char data);
extern int input_rt17  (raw_t *raw, unsigned char data);
extern int input_sbf   (raw_t *raw, unsigned char data);
extern int input_lexr  (raw_t *raw, unsigned char data);
extern int input_oem4f (raw_t *raw, FILE *fp);
extern int input_oem3f (raw_t *raw, FILE *fp);
extern int input_ubxf  (raw_t *raw, FILE *fp);
extern int input_ss2f  (raw_t *raw, FILE *fp);
extern int input_cresf (raw_t *raw, FILE *fp);
extern int input_stqf  (raw_t *raw, FILE *fp);
extern int input_gw10f (raw_t *raw, FILE *fp);
extern int input_javadf(raw_t *raw, FILE *fp);
extern int input_nvsf  (raw_t *raw, FILE *fp);
extern int input_bnxf  (raw_t *raw, FILE *fp);
extern int input_rt17f (raw_t *raw, FILE *fp);
extern int input_sbff  (raw_t *raw, FILE *fp);
extern int input_lexrf (raw_t *raw, FILE *fp);

extern int gen_ubx (const char *msg, unsigned char *buff);
extern int gen_stq (const char *msg, unsigned char *buff);
extern int gen_nvs (const char *msg, unsigned char *buff);
extern int gen_lexr(const char *msg, unsigned char *buff);


#ifdef __cplusplus
}
#endif
#endif 