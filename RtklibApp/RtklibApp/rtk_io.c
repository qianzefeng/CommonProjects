#include "rtk_io.h"
#include "rtk_trace.h"
#include "rtk_time.h"
#include "rtk_positioning.h"
#include "rtk_sat.h"

/* read station positions ------------------------------------------------------
* read positions from station position file
* args   : char  *file      I   station position file containing
*                               lat(deg) lon(deg) height(m) name in a line
*          char  *rcvs      I   station name
*          double *pos      O   station position {lat,lon,h} (rad/m)
*                               (all 0 if search error)
* return : none
*-----------------------------------------------------------------------------*/
extern void readpos(const char *file, const char *rcv, double *pos)
{
	static double poss[2048][3];
	static char stas[2048][16];
	FILE *fp;
	int i,j,len,np=0;
	char buff[256],str[256];

	trace(3,"readpos: file=%s\n",file);

	if (!(fp=fopen(file,"r"))) 
	{
		fprintf(stderr,"reference position file open error : %s\n",file);
		return;
	}
	while (np<2048&&fgets(buff,sizeof(buff),fp))
	{
		if (buff[0]=='%'||buff[0]=='#') continue;
		if (sscanf(buff,"%lf %lf %lf %s",&poss[np][0],&poss[np][1],&poss[np][2],
			str)<4) continue;
		strncpy(stas[np],str,15); stas[np++][15]='\0';
	}
	fclose(fp);
	len=(int)strlen(rcv);
	for (i=0;i<np;i++) 
	{
		if (strncmp(stas[i],rcv,len)) continue;
		for (j=0;j<3;j++) pos[j]=poss[i][j];
		pos[0]*=D2R; pos[1]*=D2R;
		return;
	}
	pos[0]=pos[1]=pos[2]=0.0;
}

/* compare observation data -------------------------------------------------*/
static int cmpobs(const void *p1, const void *p2)
{
	obsd_t *q1=(obsd_t *)p1,*q2=(obsd_t *)p2;
	double tt=timediff(q1->time,q2->time);
	if (fabs(tt)>DTTOL) return tt<0?-1:1;
	if (q1->rcv!=q2->rcv) return (int)q1->rcv-(int)q2->rcv;
	return (int)q1->sat-(int)q2->sat;
}

/* sort and unique observation data --------------------------------------------
* sort and unique observation data by time, rcv, sat
* args   : obs_t *obs    IO     observation data
* return : number of epochs
*-----------------------------------------------------------------------------*/
extern int sortobs(obs_t *obs)
{
	int i,j,n;

	trace(3,"sortobs: nobs=%d\n",obs->n);

	if (obs->n<=0) return 0;

	qsort(obs->data, obs->n, sizeof(obsd_t), cmpobs);//≈≈–Ú

	/* delete duplicated data */
	for (i=j=0;i<obs->n;i++) 
	{
		if (obs->data[i].sat!=obs->data[j].sat||
			obs->data[i].rcv!=obs->data[j].rcv||
			timediff(obs->data[i].time,obs->data[j].time)!=0.0) 
		{
			obs->data[++j]=obs->data[i];
		}
	}
	obs->n=j+1;

	for (i=n=0;i<obs->n;i=j,n++) 
	{
		for (j=i+1;j<obs->n;j++) 
		{
			if (timediff(obs->data[j].time,obs->data[i].time)>DTTOL) break;
		}
	}
	return n;
}



/* compare ephemeris ---------------------------------------------------------*/
static int cmpeph(const void *p1, const void *p2)
{
	eph_t *q1=(eph_t *)p1,*q2=(eph_t *)p2;
	return q1->ttr.time!=q2->ttr.time?(int)(q1->ttr.time-q2->ttr.time):
		(q1->toe.time!=q2->toe.time?(int)(q1->toe.time-q2->toe.time):
		q1->sat-q2->sat);
}

/* sort and unique ephemeris -------------------------------------------------*/
static void uniqeph(nav_t *nav)
{
	eph_t *nav_eph;
	int i,j;

	trace(3,"uniqeph: n=%d\n",nav->n);

	if (nav->n<=0) return;

	qsort(nav->eph,nav->n,sizeof(eph_t),cmpeph);

	for (i=1,j=0;i<nav->n;i++) 
	{
		if (nav->eph[i].sat!=nav->eph[j].sat||
			nav->eph[i].iode!=nav->eph[j].iode) 
		{
			nav->eph[++j]=nav->eph[i];
		}
	}
	nav->n=j+1;

	if (!(nav_eph=(eph_t *)realloc(nav->eph,sizeof(eph_t)*nav->n))) 
	{
		trace(1,"uniqeph malloc error n=%d\n",nav->n);
		free(nav->eph); 
		nav->eph=NULL; 
		nav->n=nav->nmax=0;
		return;
	}
	nav->eph=nav_eph;
	nav->nmax=nav->n;
	trace(4,"uniqeph: n=%d\n",nav->n);
}

/* compare glonass ephemeris -------------------------------------------------*/
static int cmpgeph(const void *p1, const void *p2)
{
	geph_t *q1=(geph_t *)p1,*q2=(geph_t *)p2;
	return q1->tof.time!=q2->tof.time?(int)(q1->tof.time-q2->tof.time):
		(q1->toe.time!=q2->toe.time?(int)(q1->toe.time-q2->toe.time):
		q1->sat-q2->sat);
}

/* sort and unique glonass ephemeris -----------------------------------------*/
static void uniqgeph(nav_t *nav)
{
	geph_t *nav_geph;
	int i,j;

	trace(3,"uniqgeph: ng=%d\n",nav->ng);

	if (nav->ng<=0) return;

	qsort(nav->geph,nav->ng,sizeof(geph_t),cmpgeph);

	for (i=j=0;i<nav->ng;i++) 
	{
		if (nav->geph[i].sat!=nav->geph[j].sat||
			nav->geph[i].toe.time!=nav->geph[j].toe.time||
			nav->geph[i].svh!=nav->geph[j].svh) 
		{
			nav->geph[++j]=nav->geph[i];
		}
	}
	nav->ng=j+1;

	if (!(nav_geph=(geph_t *)realloc(nav->geph,sizeof(geph_t)*nav->ng))) 
	{
		trace(1,"uniqgeph malloc error ng=%d\n",nav->ng);
		free(nav->geph); nav->geph=NULL; nav->ng=nav->ngmax=0;
		return;
	}
	nav->geph=nav_geph;
	nav->ngmax=nav->ng;

	trace(4,"uniqgeph: ng=%d\n",nav->ng);
}

/* compare sbas ephemeris ----------------------------------------------------*/
static int cmpseph(const void *p1, const void *p2)
{
	seph_t *q1=(seph_t *)p1,*q2=(seph_t *)p2;
	return q1->tof.time!=q2->tof.time?(int)(q1->tof.time-q2->tof.time):
		(q1->t0.time!=q2->t0.time?(int)(q1->t0.time-q2->t0.time):
		q1->sat-q2->sat);
}
/* sort and unique sbas ephemeris --------------------------------------------*/
static void uniqseph(nav_t *nav)
{
	seph_t *nav_seph;
	int i,j;

	trace(3,"uniqseph: ns=%d\n",nav->ns);

	if (nav->ns<=0) return;

	qsort(nav->seph,nav->ns,sizeof(seph_t),cmpseph);

	for (i=j=0;i<nav->ns;i++)
	{
		if (nav->seph[i].sat!=nav->seph[j].sat||
			nav->seph[i].t0.time!=nav->seph[j].t0.time) 
		{
			nav->seph[++j]=nav->seph[i];
		}
	}
	nav->ns=j+1;

	if (!(nav_seph=(seph_t *)realloc(nav->seph,sizeof(seph_t)*nav->ns)))
	{
		trace(1,"uniqseph malloc error ns=%d\n",nav->ns);
		free(nav->seph); nav->seph=NULL; nav->ns=nav->nsmax=0;
		return;
	}
	nav->seph=nav_seph;
	nav->nsmax=nav->ns;

	trace(4,"uniqseph: ns=%d\n",nav->ns);
}

/* unique ephemerides ----------------------------------------------------------
* unique ephemerides in navigation data and update carrier wave length
* args   : nav_t *nav    IO     navigation data
* return : number of epochs
*-----------------------------------------------------------------------------*/
extern void uniqnav(nav_t *nav)
{
	int i,j;

	trace(3,"uniqnav: neph=%d ngeph=%d nseph=%d\n",nav->n,nav->ng,nav->ns);

	/* unique ephemeris */
	uniqeph (nav);
	uniqgeph(nav);
	uniqseph(nav);

	/* update carrier wave length */
	for (i=0;i<MAXSAT;i++) for (j=0;j<NFREQ;j++) {
		nav->lam[i][j]=satwavelen(i+1,j,nav);
	}
}


/* screen by time --------------------------------------------------------------
* screening by time start, time end, and time interval
* args   : gtime_t time  I      time
*          gtime_t ts    I      time start (ts.time==0:no screening by ts)
*          gtime_t te    I      time end   (te.time==0:no screening by te)
*          double  tint  I      time interval (s) (0.0:no screen by tint)
* return : 1:on condition, 0:not on condition
*-----------------------------------------------------------------------------*/
extern int screent(gtime_t time, gtime_t ts, gtime_t te, double tint)
{
	return (tint<=0.0||fmod(time2gpst(time,NULL)+DTTOL,tint)<=DTTOL*2.0)&&
		(ts.time==0||timediff(time,ts)>=-DTTOL)&&
		(te.time==0||timediff(time,te)<  DTTOL);
}


/* read/save navigation data ---------------------------------------------------
* save or load navigation data
* args   : char    file  I      file path
*          nav_t   nav   O/I    navigation data
* return : status (1:ok,0:no file)
*-----------------------------------------------------------------------------*/
extern int readnav(const char *file, nav_t *nav)
{
	FILE *fp;
	eph_t eph0={0};
	char buff[4096],*p;
	int i,sat;

	trace(3,"loadnav: file=%s\n",file);

	if (!(fp=fopen(file,"r"))) return 0;

	while (fgets(buff,sizeof(buff),fp)) 
	{
		if (!strncmp(buff,"IONUTC",6)) 
		{
			for (i=0;i<8;i++) nav->ion_gps[i]=0.0;
			for (i=0;i<4;i++) nav->utc_gps[i]=0.0;
			nav->leaps=0;
			sscanf(buff,"IONUTC,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%d",
				&nav->ion_gps[0],&nav->ion_gps[1],&nav->ion_gps[2],&nav->ion_gps[3],
				&nav->ion_gps[4],&nav->ion_gps[5],&nav->ion_gps[6],&nav->ion_gps[7],
				&nav->utc_gps[0],&nav->utc_gps[1],&nav->utc_gps[2],&nav->utc_gps[3],
				&nav->leaps);
			continue;   
		}
		if ((p=strchr(buff,','))) *p='\0'; else continue;
		if (!(sat=satid2no(buff))) continue;
		nav->eph[sat-1]=eph0;
		nav->eph[sat-1].sat=sat;
		sscanf(p+1,"%d,%d,%d,%d,%ld,%ld,%ld,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,"
			"%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%d,%d",
			&nav->eph[sat-1].iode,&nav->eph[sat-1].iodc,&nav->eph[sat-1].sva ,
			&nav->eph[sat-1].svh ,&nav->eph[sat-1].toe.time,
			&nav->eph[sat-1].toc.time,&nav->eph[sat-1].ttr.time,
			&nav->eph[sat-1].A   ,&nav->eph[sat-1].e   ,&nav->eph[sat-1].i0  ,
			&nav->eph[sat-1].OMG0,&nav->eph[sat-1].omg ,&nav->eph[sat-1].M0  ,
			&nav->eph[sat-1].deln,&nav->eph[sat-1].OMGd,&nav->eph[sat-1].idot,
			&nav->eph[sat-1].crc ,&nav->eph[sat-1].crs ,&nav->eph[sat-1].cuc ,
			&nav->eph[sat-1].cus ,&nav->eph[sat-1].cic ,&nav->eph[sat-1].cis ,
			&nav->eph[sat-1].toes,&nav->eph[sat-1].fit ,&nav->eph[sat-1].f0  ,
			&nav->eph[sat-1].f1  ,&nav->eph[sat-1].f2  ,&nav->eph[sat-1].tgd[0],
			&nav->eph[sat-1].code, &nav->eph[sat-1].flag);
	}
	fclose(fp);
	return 1;
}


extern int savenav(const char *file, const nav_t *nav)
{
	FILE *fp;
	int i;
	char id[32];

	trace(3,"savenav: file=%s\n",file);

	if (!(fp=fopen(file,"w"))) return 0;

	for (i=0;i<MAXSAT;i++) 
	{
		if (nav->eph[i].ttr.time==0) continue;
		satno2id(nav->eph[i].sat,id);
		fprintf(fp,"%s,%d,%d,%d,%d,%d,%d,%d,%.14E,%.14E,%.14E,%.14E,%.14E,%.14E,"
			"%.14E,%.14E,%.14E,%.14E,%.14E,%.14E,%.14E,%.14E,%.14E,%.14E,"
			"%.14E,%.14E,%.14E,%.14E,%.14E,%d,%d\n",
			id,nav->eph[i].iode,nav->eph[i].iodc,nav->eph[i].sva ,
			nav->eph[i].svh ,(int)nav->eph[i].toe.time,
			(int)nav->eph[i].toc.time,(int)nav->eph[i].ttr.time,
			nav->eph[i].A   ,nav->eph[i].e  ,nav->eph[i].i0  ,nav->eph[i].OMG0,
			nav->eph[i].omg ,nav->eph[i].M0 ,nav->eph[i].deln,nav->eph[i].OMGd,
			nav->eph[i].idot,nav->eph[i].crc,nav->eph[i].crs ,nav->eph[i].cuc ,
			nav->eph[i].cus ,nav->eph[i].cic,nav->eph[i].cis ,nav->eph[i].toes,
			nav->eph[i].fit ,nav->eph[i].f0 ,nav->eph[i].f1  ,nav->eph[i].f2  ,
			nav->eph[i].tgd[0],nav->eph[i].code,nav->eph[i].flag);
	}
	fprintf(fp,"IONUTC,%.14E,%.14E,%.14E,%.14E,%.14E,%.14E,%.14E,%.14E,%.14E,"
		"%.14E,%.14E,%.14E,%d",
		nav->ion_gps[0],nav->ion_gps[1],nav->ion_gps[2],nav->ion_gps[3],
		nav->ion_gps[4],nav->ion_gps[5],nav->ion_gps[6],nav->ion_gps[7],
		nav->utc_gps[0],nav->utc_gps[1],nav->utc_gps[2],nav->utc_gps[3],
		nav->leaps);

	fclose(fp);
	return 1;
}

/* free observation data -------------------------------------------------------
* free memory for observation data
* args   : obs_t *obs    IO     observation data
* return : none
*-----------------------------------------------------------------------------*/
extern void freeobs(obs_t *obs)
{
	free(obs->data); obs->data=NULL; obs->n=obs->nmax=0;
}

/* free navigation data ---------------------------------------------------------
* free memory for navigation data
* args   : nav_t *nav    IO     navigation data
*          int   opt     I      option (or of followings)
*                               (0x01: gps/qzs ephmeris, 0x02: glonass ephemeris,
*                                0x04: sbas ephemeris,   0x08: precise ephemeris,
*                                0x10: precise clock     0x20: almanac,
*                                0x40: tec data)
* return : none
*-----------------------------------------------------------------------------*/
extern void freenav(nav_t *nav, int opt)
{
	if (opt&0x01) {free(nav->eph ); nav->eph =NULL; nav->n =nav->nmax =0;}
	if (opt&0x02) {free(nav->geph); nav->geph=NULL; nav->ng=nav->ngmax=0;}
	if (opt&0x04) {free(nav->seph); nav->seph=NULL; nav->ns=nav->nsmax=0;}
	if (opt&0x08) {free(nav->peph); nav->peph=NULL; nav->ne=nav->nemax=0;}
	if (opt&0x10) {free(nav->pclk); nav->pclk=NULL; nav->nc=nav->ncmax=0;}
	if (opt&0x20) {free(nav->alm ); nav->alm =NULL; nav->na=nav->namax=0;}
	if (opt&0x40) {free(nav->tec ); nav->tec =NULL; nav->nt=nav->ntmax=0;}
	if (opt&0x80) {free(nav->fcb ); nav->fcb =NULL; nav->nf=nav->nfmax=0;}
}

/* read blq record -----------------------------------------------------------*/
static int readblqrecord(FILE *fp, double *odisp)
{
	double v[11];
	char buff[256];
	int i,n=0;

	while (fgets(buff,sizeof(buff),fp)) 
	{
		if (!strncmp(buff,"$$",2)) continue;
		if (sscanf(buff,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
			v,v+1,v+2,v+3,v+4,v+5,v+6,v+7,v+8,v+9,v+10)<11) continue;
		for (i=0;i<11;i++) odisp[n+i*6]=v[i];
		if (++n==6) return 1;
	}
	return 0;
}

/* read blq ocean tide loading parameters --------------------------------------
* read blq ocean tide loading parameters
* args   : char   *file       I   BLQ ocean tide loading parameter file
*          char   *sta        I   station name
*          double *odisp      O   ocean tide loading parameters
* return : status (1:ok,0:file open error)
*-----------------------------------------------------------------------------*/
extern int readblq(const char *file, const char *sta, double *odisp)
{
	FILE *fp;
	char buff[256],staname[32]="",name[32],*p;

	/* station name to upper case */
	sscanf(sta,"%16s",staname);
	for (p=staname;(*p=(char)toupper((int)(*p)));p++) ;

	if (!(fp=fopen(file,"r"))) 
	{
		trace(2,"blq file open error: file=%s\n",file);
		return 0;
	}
	while (fgets(buff,sizeof(buff),fp)) 
	{
		if (!strncmp(buff,"$$",2)||strlen(buff)<2) continue;

		if (sscanf(buff+2,"%16s",name)<1) continue;
		for (p=name;(*p=(char)toupper((int)(*p)));p++) ;
		if (strcmp(name,staname)) continue;

		/* read blq record */
		if (readblqrecord(fp,odisp)) 
		{
			fclose(fp);
			return 1;
		}
	}
	fclose(fp);
	trace(2,"no otl parameters: sta=%s file=%s\n",sta,file);
	return 0;
}

/* read earth rotation parameters ----------------------------------------------
* read earth rotation parameters
* args   : char   *file       I   IGS ERP file (IGS ERP ver.2)
*          erp_t  *erp        O   earth rotation parameters
* return : status (1:ok,0:file open error)
*-----------------------------------------------------------------------------*/
extern int readerp(const char *file, erp_t *erp)
{
	FILE *fp;
	erpd_t *erp_data;
	double v[14]={0};
	char buff[256];

	trace(3,"readerp: file=%s\n",file);

	if (!(fp=fopen(file,"r"))) 
	{
		trace(2,"erp file open error: file=%s\n",file);
		return 0;
	}
	while (fgets(buff,sizeof(buff),fp)) 
	{
		if (sscanf(buff,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
			v,v+1,v+2,v+3,v+4,v+5,v+6,v+7,v+8,v+9,v+10,v+11,v+12,v+13)<5) {
				continue;
		}
		if (erp->n>=erp->nmax) 
		{
			erp->nmax=erp->nmax<=0?128:erp->nmax*2;
			erp_data=(erpd_t *)realloc(erp->data,sizeof(erpd_t)*erp->nmax);
			if (!erp_data) 
			{
				free(erp->data); erp->data=NULL; erp->n=erp->nmax=0;
				fclose(fp);
				return 0;
			}
			erp->data=erp_data;
		}
		erp->data[erp->n].mjd=v[0];
		erp->data[erp->n].xp=v[1]*1E-6*AS2R;
		erp->data[erp->n].yp=v[2]*1E-6*AS2R;
		erp->data[erp->n].ut1_utc=v[3]*1E-7;
		erp->data[erp->n].lod=v[4]*1E-7;
		erp->data[erp->n].xpr=v[12]*1E-6*AS2R;
		erp->data[erp->n++].ypr=v[13]*1E-6*AS2R;
	}
	fclose(fp);
	return 1;
}

/* get earth rotation parameter values -----------------------------------------
* get earth rotation parameter values
* args   : erp_t  *erp        I   earth rotation parameters
*          gtime_t time       I   time (gpst)
*          double *erpv       O   erp values {xp,yp,ut1_utc,lod} (rad,rad,s,s/d)
* return : status (1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern int geterp(const erp_t *erp, gtime_t time, double *erpv)
{
	const double ep[]={2000,1,1,12,0,0};
	double mjd,day,a;
	int i=0,j,k;

	trace(4,"geterp:\n");

	if (erp->n<=0) return 0;

	mjd=51544.5+(timediff(gpst2utc(time),epoch2time(ep)))/86400.0;

	if (mjd<=erp->data[0].mjd) 
	{
		day=mjd-erp->data[0].mjd;
		erpv[0]=erp->data[0].xp     +erp->data[0].xpr*day;
		erpv[1]=erp->data[0].yp     +erp->data[0].ypr*day;
		erpv[2]=erp->data[0].ut1_utc-erp->data[0].lod*day;
		erpv[3]=erp->data[0].lod;
		return 1;
	}
	if (mjd>=erp->data[erp->n-1].mjd) 
	{
		day=mjd-erp->data[erp->n-1].mjd;
		erpv[0]=erp->data[erp->n-1].xp     +erp->data[erp->n-1].xpr*day;
		erpv[1]=erp->data[erp->n-1].yp     +erp->data[erp->n-1].ypr*day;
		erpv[2]=erp->data[erp->n-1].ut1_utc-erp->data[erp->n-1].lod*day;
		erpv[3]=erp->data[erp->n-1].lod;
		return 1;
	}
	for (j=0,k=erp->n-1;j<k-1;) 
	{
		i=(j+k)/2;
		if (mjd<erp->data[i].mjd) k=i; else j=i;
	}
	if (erp->data[j].mjd==erp->data[j+1].mjd) 
	{
		a=0.5;
	}
	else 
	{
		a=(mjd-erp->data[j].mjd)/(erp->data[j+1].mjd-erp->data[j].mjd);
	}
	erpv[0]=(1.0-a)*erp->data[j].xp     +a*erp->data[j+1].xp;
	erpv[1]=(1.0-a)*erp->data[j].yp     +a*erp->data[j+1].yp;
	erpv[2]=(1.0-a)*erp->data[j].ut1_utc+a*erp->data[j+1].ut1_utc;
	erpv[3]=(1.0-a)*erp->data[j].lod    +a*erp->data[j+1].lod;
	return 1;
}