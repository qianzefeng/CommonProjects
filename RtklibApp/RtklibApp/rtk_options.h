#ifndef _RTKOPTION_H_
#define _RTKOPTION_H_

#include "rtk_common.h"

#ifdef __cplusplus
extern "C" {
#endif


/* options functions ---------------------------------------------------------*/
extern opt_t *searchopt(const char *name, const opt_t *opts);
extern int str2opt(opt_t *opt, const char *str);
extern int opt2str(const opt_t *opt, char *str);
extern int opt2buf(const opt_t *opt, char *buff);
extern int loadopts(const char *file, opt_t *opts);
extern int saveopts(const char *file, const char *mode, const char *comment,
					const opt_t *opts);
extern void resetsysopts(void);
extern void getsysopts(prcopt_t *popt, solopt_t *sopt, filopt_t *fopt);
extern void setsysopts(const prcopt_t *popt, const solopt_t *sopt,
					   const filopt_t *fopt);


#ifdef __cplusplus
}
#endif
#endif 