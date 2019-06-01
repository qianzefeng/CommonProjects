#ifndef _RTKSTREAM_H_
#define _RTKSTREAM_H_

#include "rtk_common.h"

#ifdef __cplusplus
extern "C" {
#endif

/* stream data input and output functions ------------------------------------*/
extern void strinitcom(void);
extern void strinit  (stream_t *stream);
extern void strlock  (stream_t *stream);
extern void strunlock(stream_t *stream);
extern int  stropen  (stream_t *stream, int type, int mode, const char *path);
extern void strclose (stream_t *stream);
extern int  strread  (stream_t *stream, unsigned char *buff, int n);
extern int  strwrite (stream_t *stream, unsigned char *buff, int n);
extern void strsync  (stream_t *stream1, stream_t *stream2);
extern int  strstat  (stream_t *stream, char *msg);
extern void strsum   (stream_t *stream, int *inb, int *inr, int *outb, int *outr);
extern void strsetopt(const int *opt);
extern gtime_t strgettime(stream_t *stream);
extern void strsendnmea(stream_t *stream, const double *pos);
extern void strsendcmd(stream_t *stream, const char *cmd);
extern void strsettimeout(stream_t *stream, int toinact, int tirecon);
extern void strsetdir(const char *dir);
extern void strsetproxy(const char *addr);


#ifdef __cplusplus
}
#endif
#endif 