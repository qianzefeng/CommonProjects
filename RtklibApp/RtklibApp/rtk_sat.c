#include "rtk_sat.h"
#include "rtk_trace.h"
#include "rtk_time.h"
#include "rtk_coord.h"
#include "rtk_sat.h"

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

static char codepris[6][MAXFREQ][16]={  /* code priority table */

	/* L1,G1E1a   L2,G2,B1     L5,G3,E5a L6,LEX,B3 E5a,B2    E5a+b */
	{"CPYWMNSL","PYWCMNDSLX","IQX"     ,""       ,""       ,""   }, /* GPS */
	{"PC"      ,"PC"        ,"IQX"     ,""       ,""       ,""   }, /* GLO */
	{"CABXZ"   ,""          ,"IQX"     ,"ABCXZ"  ,"IQX"    ,"IQX"}, /* GAL */
	{"CSLXZ"   ,"SLX"       ,"IQX"     ,"SLX"    ,""       ,""   }, /* QZS */
	{"C"       ,""          ,"IQX"     ,""       ,""       ,""   }, /* SBS */
	{"IQX"     ,"IQX"       ,"IQX"     ,"IQX"    ,"IQX"    ,""   }  /* BDS */
};

/* satellite system+prn/slot number to satellite number ------------------------
* convert satellite system+prn/slot number to satellite number
* args   : int    sys       I   satellite system (SYS_GPS,SYS_GLO,...)
*          int    prn       I   satellite prn/slot number
* return : satellite number (0:error)
*-----------------------------------------------------------------------------*/
extern int satno(int sys, int prn)
{
	if (prn<=0) return 0;
	switch (sys) 
	{
		case SYS_GPS:
			if (prn<MINPRNGPS||MAXPRNGPS<prn) return 0;
			return prn-MINPRNGPS+1;
		case SYS_GLO:
			if (prn<MINPRNGLO||MAXPRNGLO<prn) return 0;
			return NSATGPS+prn-MINPRNGLO+1;
		case SYS_GAL:
			if (prn<MINPRNGAL||MAXPRNGAL<prn) return 0;
			return NSATGPS+NSATGLO+prn-MINPRNGAL+1;
		case SYS_QZS:
			if (prn<MINPRNQZS||MAXPRNQZS<prn) return 0;
			return NSATGPS+NSATGLO+NSATGAL+prn-MINPRNQZS+1;
		case SYS_CMP:
			if (prn<MINPRNCMP||MAXPRNCMP<prn) return 0;
			return NSATGPS+NSATGLO+NSATGAL+NSATQZS+prn-MINPRNCMP+1;
		case SYS_LEO:
			if (prn<MINPRNLEO||MAXPRNLEO<prn) return 0;
			return NSATGPS+NSATGLO+NSATGAL+NSATQZS+NSATCMP+prn-MINPRNLEO+1;
		case SYS_SBS:
			if (prn<MINPRNSBS||MAXPRNSBS<prn) return 0;
			return NSATGPS+NSATGLO+NSATGAL+NSATQZS+NSATCMP+NSATLEO+prn-MINPRNSBS+1;
	}
	return 0;
}


/* satellite number to satellite system ----------------------------------------
* convert satellite number to satellite system
* args   : int    sat       I   satellite number (1-MAXSAT)
*          int    *prn      IO  satellite prn/slot number (NULL: no output)
* return : satellite system (SYS_GPS,SYS_GLO,...)
*-----------------------------------------------------------------------------*/
extern int satsys(int sat, int *prn)
{
	int sys = SYS_NONE;
	if (sat <= 0|| MAXSAT< sat) sat = 0;
	else if (sat <= NSATGPS) 
	{
		sys = SYS_GPS; 
		sat += MINPRNGPS-1;
	}
	else if ((sat -= NSATGPS)<= NSATGLO) 
	{
		sys = SYS_GLO; 
		sat += MINPRNGLO-1;
	}
	else if ((sat-=NSATGLO)<= NSATGAL) 
	{
		sys = SYS_GAL; 
		sat += MINPRNGAL-1;
	}
	else if ((sat-=NSATGAL)<=NSATQZS) 
	{
		sys = SYS_QZS; 
		sat += MINPRNQZS-1; 
	}
	else if ((sat-=NSATQZS)<=NSATCMP)
	{
		sys = SYS_CMP; 
		sat += MINPRNCMP-1;
	}
	else if ((sat -= NSATCMP) <= NSATLEO) 
	{
		sys = SYS_LEO; 
		sat += MINPRNLEO-1; 
	}
	else if ((sat-=NSATLEO) <= NSATSBS) 
	{
		sys = SYS_SBS; 
		sat += MINPRNSBS-1; 
	}
	else sat=0;
	if (prn) *prn=sat;
	return sys;
}

/* satellite id to satellite number --------------------------------------------
* convert satellite id to satellite number
* args   : char   *id       I   satellite id (nn,Gnn,Rnn,Enn,Jnn,Cnn or Snn)
* return : satellite number (0: error)
* notes  : 120-138 and 193-195 are also recognized as sbas and qzss
*-----------------------------------------------------------------------------*/
extern int satid2no(const char *id)
{
	int sys,prn;
	char code;

	if (sscanf(id,"%d", &prn)==1) 
	{
		if      (MINPRNGPS <= prn && prn <= MAXPRNGPS) sys = SYS_GPS;
		else if (MINPRNSBS <= prn && prn <= MAXPRNSBS) sys = SYS_SBS;
		else if (MINPRNQZS <= prn && prn <= MAXPRNQZS) sys = SYS_QZS;
		else return 0;
		return satno(sys,prn);
	}
	if (sscanf(id, "%c%d", &code, &prn)<2) return 0;

	switch (code)
	{
		case 'G': sys=SYS_GPS; prn+=MINPRNGPS-1; break;
		case 'R': sys=SYS_GLO; prn+=MINPRNGLO-1; break;
		case 'E': sys=SYS_GAL; prn+=MINPRNGAL-1; break;
		case 'J': sys=SYS_QZS; prn+=MINPRNQZS-1; break;
		case 'C': sys=SYS_CMP; prn+=MINPRNCMP-1; break;
		case 'L': sys=SYS_LEO; prn+=MINPRNLEO-1; break;
		case 'S': sys=SYS_SBS; prn+=100; break;
		default: return 0;
	}
	return satno(sys, prn);
}

/* satellite number to satellite id --------------------------------------------
* convert satellite number to satellite id
* args   : int    sat       I   satellite number
*          char   *id       O   satellite id (Gnn,Rnn,Enn,Jnn,Cnn or nnn)
* return : none
*-----------------------------------------------------------------------------*/
extern void satno2id(int sat, char *id)
{
	int prn;
	switch (satsys(sat,&prn)) 
	{
		case SYS_GPS: sprintf(id,"G%02d",prn-MINPRNGPS+1); return;
		case SYS_GLO: sprintf(id,"R%02d",prn-MINPRNGLO+1); return;
		case SYS_GAL: sprintf(id,"E%02d",prn-MINPRNGAL+1); return;
		case SYS_QZS: sprintf(id,"J%02d",prn-MINPRNQZS+1); return;
		case SYS_CMP: sprintf(id,"C%02d",prn-MINPRNCMP+1); return;
		case SYS_LEO: sprintf(id,"L%02d",prn-MINPRNLEO+1); return;
		case SYS_SBS: sprintf(id,"%03d" ,prn); return;
	}
	strcpy(id,"");
}


/* obs type string to obs code -------------------------------------------------
* convert obs code type string to obs code
* args   : char   *str   I      obs code string ("1C","1P","1Y",...)
*          int    *freq  IO     frequency (1:L1,2:L2,3:L5,4:L6,5:L7,6:L8,0:err)
*                               (NULL: no output)
* return : obs code (CODE_???)
* notes  : obs codes are based on reference [6] and qzss extension
*-----------------------------------------------------------------------------*/
extern unsigned char obs2code(const char *obs, int *freq)
{
	int i;
	if (freq) *freq=0;
	for (i=1;*obscodes[i];i++) 
	{
		if (strcmp(obscodes[i],obs)) continue;
		if (freq) *freq=obsfreqs[i];
		return (unsigned char)i;
	}
	return CODE_NONE;
}


/* obs code to obs code string -------------------------------------------------
* convert obs code to obs code string
* args   : unsigned char code I obs code (CODE_???)
*          int    *freq  IO     frequency (1:L1,2:L2,3:L5,4:L6,5:L7,6:L8,0:err)
*                               (NULL: no output)
* return : obs code string ("1C","1P","1P",...)
* notes  : obs codes are based on reference [6] and qzss extension
*-----------------------------------------------------------------------------*/
extern char *code2obs(unsigned char code, int *freq)
{
	if (freq) *freq=0;
	if (code<=CODE_NONE||MAXCODE<code) return "";
	if (freq) *freq=obsfreqs[code];
	return obscodes[code];
}

/* set code priority -----------------------------------------------------------
* set code priority for multiple codes in a frequency
* args   : int    sys     I     system (or of SYS_???)
*          int    freq    I     frequency (1:L1,2:L2,3:L5,4:L6,5:L7,6:L8)
*          char   *pri    I     priority of codes (series of code characters)
*                               (higher priority precedes lower)
* return : none
*-----------------------------------------------------------------------------*/
extern void setcodepri(int sys, int freq, const char *pri)
{
	trace(3,"setcodepri:sys=%d freq=%d pri=%s\n",sys,freq,pri);

	if (freq<=0||MAXFREQ<freq) return;
	if (sys&SYS_GPS) strcpy(codepris[0][freq-1],pri);
	if (sys&SYS_GLO) strcpy(codepris[1][freq-1],pri);
	if (sys&SYS_GAL) strcpy(codepris[2][freq-1],pri);
	if (sys&SYS_QZS) strcpy(codepris[3][freq-1],pri);
	if (sys&SYS_SBS) strcpy(codepris[4][freq-1],pri);
	if (sys&SYS_CMP) strcpy(codepris[5][freq-1],pri);
}


/* get code priority -----------------------------------------------------------
* get code priority for multiple codes in a frequency
* args   : int    sys     I     system (SYS_???)
*          unsigned char code I obs code (CODE_???)
*          char   *opt    I     code options (NULL:no option)
* return : priority (15:highest-1:lowest,0:error)
*-----------------------------------------------------------------------------*/
extern int getcodepri(int sys, unsigned char code, const char *opt)
{
	const char *p,*optstr;
	char *obs,str[8]="";
	int i,j;

	switch (sys) 
	{
		case SYS_GPS: i=0; optstr="-GL%2s"; break;
		case SYS_GLO: i=1; optstr="-RL%2s"; break;
		case SYS_GAL: i=2; optstr="-EL%2s"; break;
		case SYS_QZS: i=3; optstr="-JL%2s"; break;
		case SYS_SBS: i=4; optstr="-SL%2s"; break;
		case SYS_CMP: i=5; optstr="-CL%2s"; break;
		default: return 0;
	}
	obs=code2obs(code,&j);

	/* parse code options */
	for (p=opt;p&&(p=strchr(p,'-'));p++) 
	{
		if (sscanf(p,optstr,str)<1||str[0]!=obs[0]) continue;
		return str[1]==obs[1]?15:0;
	}
	/* search code priority */
	return (p=strchr(codepris[i][j-1],obs[1]))?14-(int)(p-codepris[i][j-1]):0;
}


/* test excluded satellite -----------------------------------------------------
* test excluded satellite
* args   : int    sat       I   satellite number
*          int    svh       I   sv health flag
*          prcopt_t *opt    I   processing options (NULL: not used)
* return : status (1:excluded,0:not excluded)
*-----------------------------------------------------------------------------*/
extern int satexclude(int sat, int svh, const prcopt_t *opt)
{
	int sys=satsys(sat,NULL);

	if (svh<0) return 1; /* ephemeris unavailable */

	if (opt) 
	{
		if (opt->exsats[sat-1]==1) return 1; /* excluded satellite */
		if (opt->exsats[sat-1]==2) return 0; /* included satellite */
		if (!(sys&opt->navsys)) return 1; /* unselected sat sys */
	}
	if (sys==SYS_QZS) svh&=0xFE; /* mask QZSS LEX health */
	if (svh) 
	{
		trace(3,"unhealthy satellite: sat=%3d svh=%02X\n",sat,svh);
		return 1;
	}
	return 0;
}


/* test SNR mask ---------------------------------------------------------------
* test SNR mask
* args   : int    base      I   rover or base-station (0:rover,1:base station)
*          int    freq      I   frequency (0:L1,1:L2,2:L3,...)
*          double el        I   elevation angle (rad)
*          double snr       I   C/N0 (dBHz)
*          snrmask_t *mask  I   SNR mask
* return : status (1:masked,0:unmasked)
*-----------------------------------------------------------------------------*/
extern int testsnr(int base, int freq, double el, double snr,
				   const snrmask_t *mask)
{
	double minsnr,a;
	int i;

	if (!mask->ena[base]||freq<0||freq>=NFREQ) return 0;

	a=(el*R2D+5.0)/10.0;
	i=(int)floor(a); a-=i;
	if      (i<1) minsnr=mask->mask[freq][0];
	else if (i>8) minsnr=mask->mask[freq][8];
	else minsnr=(1.0-a)*mask->mask[freq][i-1]+a*mask->mask[freq][i];

	return snr<minsnr;
}