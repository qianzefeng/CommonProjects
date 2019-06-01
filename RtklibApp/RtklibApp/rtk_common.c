#include "rtk_common.h"


const double chisqr[100]={      /* chi-sqr(n) (alpha=0.001) */
	10.8,13.8,16.3,18.5,20.5,22.5,24.3,26.1,27.9,29.6,
	31.3,32.9,34.5,36.1,37.7,39.3,40.8,42.3,43.8,45.3,
	46.8,48.3,49.7,51.2,52.6,54.1,55.5,56.9,58.3,59.7,
	61.1,62.5,63.9,65.2,66.6,68.0,69.3,70.7,72.1,73.4,
	74.7,76.0,77.3,78.6,80.0,81.3,82.6,84.0,85.4,86.7,
	88.0,89.3,90.6,91.9,93.3,94.7,96.0,97.4,98.7,100 ,
	101 ,102 ,103 ,104 ,105 ,107 ,108 ,109 ,110 ,112 ,
	113 ,114 ,115 ,116 ,118 ,119 ,120 ,122 ,123 ,125 ,
	126 ,127 ,128 ,129 ,131 ,132 ,133 ,134 ,135 ,137 ,
	138 ,139 ,140 ,142 ,143 ,144 ,145 ,147 ,148 ,149
};

const double lam_carr[]={       /* carrier wave length (m) */
	CLIGHT/FREQ1,CLIGHT/FREQ2,CLIGHT/FREQ5,CLIGHT/FREQ6,CLIGHT/FREQ7,CLIGHT/FREQ8
};

const prcopt_t prcopt_default={ /* defaults processing options */
	PMODE_SINGLE,0,2,SYS_GPS,   /* mode,soltype,nf,navsys */
	15.0*D2R,{{0,0}},           /* elmin,snrmask */
	0,1,1,1,                    /* sateph,modear,glomodear,bdsmodear */
	5,0,10,1,                   /* maxout,minlock,minfix,armaxiter */
	0,0,0,0,                    /* estion,esttrop,dynamics,tidecorr */
	1,0,0,0,0,                  /* niter,codesmooth,intpref,sbascorr,sbassatsel */
	0,0,                        /* rovpos,refpos */
	{100.0,100.0},              /* eratio[] */
	{100.0,0.003,0.003,0.0,1.0}, /* err[] */
	{30.0,0.03,0.3},            /* std[] */
	{1E-4,1E-3,1E-4,1E-1,1E-2,0.0}, /* prn[] */
	5E-12,                      /* sclkstab */
	{3.0,0.9999,0.25,0.1,0.05}, /* thresar */
	0.0,0.0,0.05,               /* elmaskar,almaskhold,thresslip */
	30.0,30.0,30.0,             /* maxtdif,maxinno,maxgdop */
	{0},{0},{0},                /* baseline,ru,rb */
	{"",""},                    /* anttype */
	{{0}},{{0}},{0}             /* antdel,pcv,exsats */
};
const solopt_t solopt_default={ /* defaults solution output options */
	SOLF_LLH,TIMES_GPST,1,3,    /* posf,times,timef,timeu */
	0,1,0,0,0,0,                /* degf,outhead,outopt,datum,height,geoid */
	0,0,0,                      /* solstatic,sstat,trace */
	{0.0,0.0},                  /* nmeaintv */
	" ",""                      /* separator/program name */
};
const char *formatstrs[]={      /* stream format strings */
	"RTCM 2",                   /*  0 */
	"RTCM 3",                   /*  1 */
	"NovAtel OEM6",             /*  2 */
	"NovAtel OEM3",             /*  3 */
	"u-blox",                   /*  4 */
	"Superstar II",             /*  5 */
	"Hemisphere",               /*  6 */
	"SkyTraq",                  /*  7 */
	"GW10",                     /*  8 */
	"Javad",                    /*  9 */
	"NVS BINR",                 /* 10 */
	"BINEX",                    /* 11 */
	"Trimble RT17",             /* 12 */
	"Septentrio",               /* 13 */
	"LEX Receiver",             /* 14 */
	"RINEX",                    /* 15 */
	"SP3",                      /* 16 */
	"RINEX CLK",                /* 17 */
	"SBAS",                     /* 18 */
	"NMEA 0183",                /* 19 */
	NULL
};
static char *obscodes[]={       /* observation code strings */

	""  ,"1C","1P","1W","1Y", "1M","1N","1S","1L","1E", /*  0- 9 */
	"1A","1B","1X","1Z","2C", "2D","2S","2L","2X","2P", /* 10-19 */
	"2W","2Y","2M","2N","5I", "5Q","5X","7I","7Q","7X", /* 20-29 */
	"6A","6B","6C","6X","6Z", "6S","6L","8L","8Q","8X", /* 30-39 */
	"2I","2Q","6I","6Q","3I", "3Q","3X","1I","1Q",""    /* 40-49 */
};
static unsigned char obsfreqs[]={ /* 1:L1,2:L2,3:L5,4:L6,5:L7,6:L8,7:L3 */

	0, 1, 1, 1, 1,  1, 1, 1, 1, 1, /*  0- 9 */
	1, 1, 1, 1, 2,  2, 2, 2, 2, 2, /* 10-19 */
	2, 2, 2, 2, 3,  3, 3, 5, 5, 5, /* 20-29 */
	4, 4, 4, 4, 4,  4, 4, 6, 6, 6, /* 30-39 */
	2, 2, 4, 4, 3,  3, 3, 1, 1, 0  /* 40-49 */
};