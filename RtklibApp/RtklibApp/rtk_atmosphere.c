#include "rtk_atmosphere.h"
#include "rtk_mat.h"
#include "rtk_time.h"
#include "rtk_trace.h"
#include "rtk_platform.h"
#include "rtk_sat.h"
#include "rtk_sbas.h"
#include "rtk_qzss.h"

#define SQR(x)      ((x)*(x))
#define VAR_NOTEC   SQR(30.0)   /* variance of no tec */
#define MIN_EL      0.0         /* min elevation angle (rad) */
#define MIN_HGT     -1000.0     /* min user height (m) */

#define MAXITR      10          /* max number of iteration for point pos */
#define ERR_ION     5.0         /* ionospheric delay std (m) */
#define ERR_TROP    3.0         /* tropspheric delay std (m) */
#define ERR_SAAS    0.3         /* saastamoinen model error std (m) */
#define ERR_BRDCI   0.5         /* broadcast iono model error factor */
#define ERR_CBIAS   0.3         /* code bias error std (m) */
#define REL_HUMI    0.7         /* relative humidity for saastamoinen model */

/* get index -----------------------------------------------------------------*/
static int getindex(double value, const double *range)
{
	if (range[2]==0.0) return 0;
	if (range[1]>0.0&&(value<range[0]||range[1]<value)) return -1;
	if (range[1]<0.0&&(value<range[1]||range[0]<value)) return -1;
	return (int)floor((value-range[0])/range[2]+0.5);
}

/* get number of items -------------------------------------------------------*/
static int nitem(const double *range)
{
	return getindex(range[1],range)+1;
}

/* data index (i:lat,j:lon,k:hgt) --------------------------------------------*/
static int dataindex(int i, int j, int k, const int *ndata)
{
	if (i<0||ndata[0]<=i||j<0||ndata[1]<=j||k<0||ndata[2]<=k) return -1;
	return i+ndata[0]*(j+ndata[1]*k);
}

/* add tec data to navigation data -------------------------------------------*/
static tec_t *addtec(const double *lats, const double *lons, const double *hgts, double rb, nav_t *nav)
{
	tec_t *p,*nav_tec;
	gtime_t time0={0};
	int i,n,ndata[3];

	trace(3,"addtec  :\n");

	ndata[0]=nitem(lats);
	ndata[1]=nitem(lons);
	ndata[2]=nitem(hgts);
	if (ndata[0]<=1||ndata[1]<=1||ndata[2]<=0) return NULL;

	if (nav->nt>=nav->ntmax) 
	{
		nav->ntmax+=256;
		if (!(nav_tec=(tec_t *)realloc(nav->tec,sizeof(tec_t)*nav->ntmax))) 
		{
			trace(1,"readionex malloc error ntmax=%d\n",nav->ntmax);
			free(nav->tec); nav->tec=NULL; nav->nt=nav->ntmax=0;
			return NULL;
		}
		nav->tec=nav_tec;
	}
	p=nav->tec+nav->nt;
	p->time=time0;
	p->rb=rb;
	for (i=0;i<3;i++) 
	{
		p->ndata[i]=ndata[i];
		p->lats[i]=lats[i];
		p->lons[i]=lons[i];
		p->hgts[i]=hgts[i];
	}
	n=ndata[0]*ndata[1]*ndata[2];

	if (!(p->data=(double *)malloc(sizeof(double)*n))||
		!(p->rms =(float  *)malloc(sizeof(float )*n))) 
	{
			return NULL;
	}
	for (i=0;i<n;i++) 
	{
		p->data[i]=0.0;
		p->rms [i]=0.0f;
	}
	nav->nt++;
	return p;
}


/* ionosphere model ------------------------------------------------------------
* compute ionospheric delay by broadcast ionosphere model (klobuchar model)
* args   : gtime_t t        I   time (gpst)
*          double *ion      I   iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3}
*          double *pos      I   receiver position {lat,lon,h} (rad,m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
* return : ionospheric delay (L1) (m)
*-----------------------------------------------------------------------------*/
extern double ionmodel(gtime_t t, const double *ion, const double *pos,
					   const double *azel)
{
	const double ion_default[]={ /* 2004/1/1 */
		0.1118E-07,-0.7451E-08,-0.5961E-07, 0.1192E-06,
		0.1167E+06,-0.2294E+06,-0.1311E+06, 0.1049E+07
	};
	double tt,f,psi,phi,lam,amp,per,x;
	int week;

	if (pos[2]<-1E3||azel[1]<=0) return 0.0;
	if (norm(ion,8)<=0.0) ion=ion_default;

	/* earth centered angle (semi-circle) */
	psi=0.0137/(azel[1]/PI+0.11)-0.022;

	/* subionospheric latitude/longitude (semi-circle) */
	phi=pos[0]/PI+psi*cos(azel[0]);
	if      (phi> 0.416) phi= 0.416;
	else if (phi<-0.416) phi=-0.416;
	lam=pos[1]/PI+psi*sin(azel[0])/cos(phi*PI);

	/* geomagnetic latitude (semi-circle) */
	phi+=0.064*cos((lam-1.617)*PI);

	/* local time (s) */
	tt=43200.0*lam+time2gpst(t,&week);
	tt-=floor(tt/86400.0)*86400.0; /* 0<=tt<86400 */

	/* slant factor */
	f=1.0+16.0*pow(0.53-azel[1]/PI,3.0);

	/* ionospheric delay */
	amp=ion[0]+phi*(ion[1]+phi*(ion[2]+phi*ion[3]));
	per=ion[4]+phi*(ion[5]+phi*(ion[6]+phi*ion[7]));
	amp=amp<    0.0?    0.0:amp;
	per=per<72000.0?72000.0:per;
	x=2.0*PI*(tt-50400.0)/per;

	return CLIGHT*f*(fabs(x)<1.57?5E-9+amp*(1.0+x*x*(-0.5+x*x/24.0)):5E-9);
}

/* ionosphere mapping function -------------------------------------------------
* compute ionospheric delay mapping function by single layer model
* args   : double *pos      I   receiver position {lat,lon,h} (rad,m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
* return : ionospheric mapping function
*-----------------------------------------------------------------------------*/
extern double ionmapf(const double *pos, const double *azel)
{
	if (pos[2]>=HION) return 1.0;
	return 1.0/cos(asin((RE_WGS84+pos[2])/(RE_WGS84+HION)*sin(PI/2.0-azel[1])));
}


/* ionospheric pierce point position -------------------------------------------
* compute ionospheric pierce point (ipp) position and slant factor
* args   : double *pos      I   receiver position {lat,lon,h} (rad,m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
*          double re        I   earth radius (km)
*          double hion      I   altitude of ionosphere (km)
*          double *posp     O   pierce point position {lat,lon,h} (rad,m)
* return : slant factor
* notes  : see ref [2], only valid on the earth surface
*          fixing bug on ref [2] A.4.4.10.1 A-22,23
*-----------------------------------------------------------------------------*/
extern double ionppp(const double *pos, const double *azel, double re,
					 double hion, double *posp)
{
	double cosaz,rp,ap,sinap,tanap;

	rp=re/(re+hion)*cos(azel[1]);
	ap=PI/2.0-azel[1]-asin(rp);
	sinap=sin(ap);
	tanap=tan(ap);
	cosaz=cos(azel[0]);
	posp[0]=asin(sin(pos[0])*cos(ap)+cos(pos[0])*sinap*cosaz);

	if ((pos[0]> 70.0*D2R&& tanap*cosaz>tan(PI/2.0-pos[0]))||
		(pos[0]<-70.0*D2R&&-tanap*cosaz>tan(PI/2.0+pos[0]))) {
			posp[1]=pos[1]+PI-asin(sinap*sin(azel[0])/cos(posp[0]));
	}
	else {
		posp[1]=pos[1]+asin(sinap*sin(azel[0])/cos(posp[0]));
	}
	return 1.0/sqrt(1.0-rp*rp);
}


/* troposphere model -----------------------------------------------------------
* compute tropospheric delay by standard atmosphere and saastamoinen model
* args   : gtime_t time     I   time
*          double *pos      I   receiver position {lat,lon,h} (rad,m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
*          double humi      I   relative humidity
* return : tropospheric delay (m)
*-----------------------------------------------------------------------------*/
extern double tropmodel(gtime_t time, const double *pos, const double *azel,
						double humi)
{
	const double temp0=15.0; /* temparature at sea level */
	double hgt,pres,temp,e,z,trph,trpw;

	if (pos[2]<-100.0||1E4<pos[2]||azel[1]<=0) return 0.0;

	/* standard atmosphere */
	hgt=pos[2]<0.0?0.0:pos[2];

	pres=1013.25*pow(1.0-2.2557E-5*hgt,5.2568);
	temp=temp0-6.5E-3*hgt+273.16;
	e=6.108*humi*exp((17.15*temp-4684.0)/(temp-38.45));

	/* saastamoninen model */
	z=PI/2.0-azel[1];
	trph=0.0022768*pres/(1.0-0.00266*cos(2.0*pos[0])-0.00028*hgt/1E3)/cos(z);
	trpw=0.002277*(1255.0/temp+0.05)*e/cos(z);
	return trph+trpw;
}

static double interpc(const double coef[], double lat)
{
	int i=(int)(lat/15.0);
	if (i<1) return coef[0]; else if (i>4) return coef[4];
	return coef[i-1]*(1.0-lat/15.0+i)+coef[i]*(lat/15.0-i);
}

static double mapf(double el, double a, double b, double c)
{
	double sinel=sin(el);
	return (1.0+a/(1.0+b/(1.0+c)))/(sinel+(a/(sinel+b/(sinel+c))));
}

static double nmf(gtime_t time, const double pos[], const double azel[], double *mapfw)
{
	/* ref [5] table 3 */
	/* hydro-ave-a,b,c, hydro-amp-a,b,c, wet-a,b,c at latitude 15,30,45,60,75 */
	const double coef[][5]={
		{ 1.2769934E-3, 1.2683230E-3, 1.2465397E-3, 1.2196049E-3, 1.2045996E-3},
		{ 2.9153695E-3, 2.9152299E-3, 2.9288445E-3, 2.9022565E-3, 2.9024912E-3},
		{ 62.610505E-3, 62.837393E-3, 63.721774E-3, 63.824265E-3, 64.258455E-3},

		{ 0.0000000E-0, 1.2709626E-5, 2.6523662E-5, 3.4000452E-5, 4.1202191E-5},
		{ 0.0000000E-0, 2.1414979E-5, 3.0160779E-5, 7.2562722E-5, 11.723375E-5},
		{ 0.0000000E-0, 9.0128400E-5, 4.3497037E-5, 84.795348E-5, 170.37206E-5},

		{ 5.8021897E-4, 5.6794847E-4, 5.8118019E-4, 5.9727542E-4, 6.1641693E-4},
		{ 1.4275268E-3, 1.5138625E-3, 1.4572752E-3, 1.5007428E-3, 1.7599082E-3},
		{ 4.3472961E-2, 4.6729510E-2, 4.3908931E-2, 4.4626982E-2, 5.4736038E-2}
	};
	const double aht[]={ 2.53E-5, 5.49E-3, 1.14E-3}; /* height correction */

	double y,cosy,ah[3],aw[3],dm,el=azel[1],lat=pos[0]*R2D,hgt=pos[2];
	int i;

	if (el<=0.0) {
		if (mapfw) *mapfw=0.0;
		return 0.0;
	}
	/* year from doy 28, added half a year for southern latitudes */
	y=(time2doy(time)-28.0)/365.25+(lat<0.0?0.5:0.0);

	cosy=cos(2.0*PI*y);
	lat=fabs(lat);

	for (i=0;i<3;i++) 
	{
		ah[i]=interpc(coef[i  ],lat)-interpc(coef[i+3],lat)*cosy;
		aw[i]=interpc(coef[i+6],lat);
	}
	/* ellipsoidal height is used instead of height above sea level */
	dm=(1.0/sin(el)-mapf(el,aht[0],aht[1],aht[2]))*hgt/1E3;

	if (mapfw) *mapfw=mapf(el,aw[0],aw[1],aw[2]);

	return mapf(el,ah[0],ah[1],ah[2])+dm;
}

/* troposphere mapping function ------------------------------------------------
* compute tropospheric mapping function by NMF
* args   : gtime_t t        I   time
*          double *pos      I   receiver position {lat,lon,h} (rad,m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
*          double *mapfw    IO  wet mapping function (NULL: not output)
* return : dry mapping function
* note   : see ref [5] (NMF) and [9] (GMF)
*          original JGR paper of [5] has bugs in eq.(4) and (5). the corrected
*          paper is obtained from:
*          ftp://web.haystack.edu/pub/aen/nmf/NMF_JGR.pdf
*-----------------------------------------------------------------------------*/
extern double tropmapf(gtime_t time, const double pos[], const double azel[],
					   double *mapfw)
{
#ifdef IERS_MODEL
	const double ep[]={2000,1,1,12,0,0};
	double mjd,lat,lon,hgt,zd,gmfh,gmfw;
#endif
	trace(4,"tropmapf: pos=%10.6f %11.6f %6.1f azel=%5.1f %4.1f\n",
		pos[0]*R2D,pos[1]*R2D,pos[2],azel[0]*R2D,azel[1]*R2D);

	if (pos[2]<-1000.0||pos[2]>20000.0) {
		if (mapfw) *mapfw=0.0;
		return 0.0;
	}
#ifdef IERS_MODEL
	mjd=51544.5+(timediff(time,epoch2time(ep)))/86400.0;
	lat=pos[0];
	lon=pos[1];
	hgt=pos[2]-geoidh(pos); /* height in m (mean sea level) */
	zd =PI/2.0-azel[1];

	/* call GMF */
	gmf_(&mjd,&lat,&lon,&hgt,&zd,&gmfh,&gmfw);

	if (mapfw) *mapfw=gmfw;
	return gmfh;
#else
	return nmf(time,pos,azel,mapfw); /* NMF */
#endif
}

/* interpolate tec grid data -------------------------------------------------*/
static int interptec(const tec_t *tec, int k, const double *posp, double *value, double *rms)
{
	double dlat,dlon,a,b,d[4]={0},r[4]={0};
	int i,j,n,index;

	trace(3,"interptec: k=%d posp=%.2f %.2f\n",k,posp[0]*R2D,posp[1]*R2D);
	*value=*rms=0.0;

	if (tec->lats[2]==0.0||tec->lons[2]==0.0) return 0;

	dlat=posp[0]*R2D-tec->lats[0];
	dlon=posp[1]*R2D-tec->lons[0];
	if (tec->lons[2]>0.0) dlon-=floor( dlon/360)*360.0; /*  0<=dlon<360 */
	else                  dlon+=floor(-dlon/360)*360.0; /* -360<dlon<=0 */

	a=dlat/tec->lats[2];
	b=dlon/tec->lons[2];
	i=(int)floor(a); a-=i;
	j=(int)floor(b); b-=j;

	/* get gridded tec data */
	for (n=0;n<4;n++) 
	{
		if ((index=dataindex(i+(n%2),j+(n<2?0:1),k,tec->ndata))<0) continue;
		d[n]=tec->data[index];
		r[n]=tec->rms [index];
	}
	if (d[0]>0.0&&d[1]>0.0&&d[2]>0.0&&d[3]>0.0) 
	{
		/* bilinear interpolation (inside of grid) */
		*value=(1.0-a)*(1.0-b)*d[0]+a*(1.0-b)*d[1]+(1.0-a)*b*d[2]+a*b*d[3];
		*rms  =(1.0-a)*(1.0-b)*r[0]+a*(1.0-b)*r[1]+(1.0-a)*b*r[2]+a*b*r[3];
	}
	/* nearest-neighbour extrapolation (outside of grid) */
	else if (a<=0.5&&b<=0.5&&d[0]>0.0) {*value=d[0]; *rms=r[0];}
	else if (a> 0.5&&b<=0.5&&d[1]>0.0) {*value=d[1]; *rms=r[1];}
	else if (a<=0.5&&b> 0.5&&d[2]>0.0) {*value=d[2]; *rms=r[2];}
	else if (a> 0.5&&b> 0.5&&d[3]>0.0) {*value=d[3]; *rms=r[3];}
	else {
		i=0;
		for (n=0;n<4;n++) if (d[n]>0.0) {i++; *value+=d[n]; *rms+=r[n];}
		if(i==0) return 0;
		*value/=i; *rms/=i;
	}
	return 1;
}

/* ionosphere delay by tec grid data -----------------------------------------*/
static int iondelay(gtime_t time, const tec_t *tec, const double *pos,
					const double *azel, int opt, double *delay, double *var)
{
	const double fact=40.30E16/FREQ1/FREQ1; /* tecu->L1 iono (m) */
	double fs,posp[3]={0},vtec,rms,hion,rp;
	int i;

	trace(3,"iondelay: time=%s pos=%.1f %.1f azel=%.1f %.1f\n",time_str(time,0),
		pos[0]*R2D,pos[1]*R2D,azel[0]*R2D,azel[1]*R2D);

	*delay=*var=0.0;

	for (i=0;i<tec->ndata[2];i++) { /* for a layer */

		hion=tec->hgts[0]+tec->hgts[2]*i;

		/* ionospheric pierce point position */
		fs=ionppp(pos,azel,tec->rb,hion,posp);

		if (opt&2) {
			/* modified single layer mapping function (M-SLM) ref [2] */
			rp=tec->rb/(tec->rb+hion)*sin(0.9782*(PI/2.0-azel[1]));
			fs=1.0/sqrt(1.0-rp*rp);
		}
		if (opt&1) {
			/* earth rotation correction (sun-fixed coordinate) */
			posp[1]+=2.0*PI*timediff(time,tec->time)/86400.0;
		}
		/* interpolate tec grid data */
		if (!interptec(tec,i,posp,&vtec,&rms)) return 0;

		*delay+=fact*fs*vtec;
		*var+=fact*fact*fs*fs*rms*rms;
	}
	trace(4,"iondelay: delay=%7.2f std=%6.2f\n",*delay,sqrt(*var));

	return 1;
}


/* ionosphere model by tec grid data -------------------------------------------
* compute ionospheric delay by tec grid data
* args   : gtime_t time     I   time (gpst)
*          nav_t  *nav      I   navigation data
*          double *pos      I   receiver position {lat,lon,h} (rad,m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
*          int    opt       I   model option
*                                bit0: 0:earth-fixed,1:sun-fixed
*                                bit1: 0:single-layer,1:modified single-layer
*          double *delay    O   ionospheric delay (L1) (m)
*          double *var      O   ionospheric dealy (L1) variance (m^2)
* return : status (1:ok,0:error)
* notes  : before calling the function, read tec grid data by calling readtec()
*          return ok with delay=0 and var=VAR_NOTEC if el<MIN_EL or h<MIN_HGT
*-----------------------------------------------------------------------------*/
extern int iontec(gtime_t time, const nav_t *nav, const double *pos,
				  const double *azel, int opt, double *delay, double *var)
{
	double dels[2],vars[2],a,tt;
	int i,stat[2];

	trace(3,"iontec  : time=%s pos=%.1f %.1f azel=%.1f %.1f\n",time_str(time,0),
		pos[0]*R2D,pos[1]*R2D,azel[0]*R2D,azel[1]*R2D);

	if (azel[1]<MIN_EL||pos[2]<MIN_HGT) 
	{
		*delay=0.0;
		*var=VAR_NOTEC;
		return 1;
	}
	for (i=0;i<nav->nt;i++) 
	{
		if (timediff(nav->tec[i].time,time)>0.0) break;
	}
	if (i==0||i>=nav->nt) 
	{
		trace(2,"%s: tec grid out of period\n",time_str(time,0));
		return 0;
	}
	if ((tt=timediff(nav->tec[i].time,nav->tec[i-1].time))==0.0) 
	{
		trace(2,"tec grid time interval error\n");
		return 0;
	}
	/* ionospheric delay by tec grid data */
	stat[0]=iondelay(time,nav->tec+i-1,pos,azel,opt,dels  ,vars  );
	stat[1]=iondelay(time,nav->tec+i  ,pos,azel,opt,dels+1,vars+1);

	if (!stat[0]&&!stat[1]) 
	{
		trace(2,"%s: tec grid out of area pos=%6.2f %7.2f azel=%6.1f %5.1f\n",
			time_str(time,0),pos[0]*R2D,pos[1]*R2D,azel[0]*R2D,azel[1]*R2D);
		return 0;
	}
	if (stat[0]&&stat[1]) 
	{ /* linear interpolation by time */
		a=timediff(time,nav->tec[i-1].time)/tt;
		*delay=dels[0]*(1.0-a)+dels[1]*a;
		*var  =vars[0]*(1.0-a)+vars[1]*a;
	}
	else if (stat[0]) 
	{ 
		/* nearest-neighbour extrapolation by time */
		*delay=dels[0];
		*var  =vars[0];
	}
	else 
	{
		*delay=dels[1];
		*var  =vars[1];
	}
	trace(3,"iontec  : delay=%5.2f std=%5.2f\n",*delay,sqrt(*var));
	return 1;
}


/* read ionex dcb aux data ----------------------------------------------------*/
static void readionexdcb(FILE *fp, double *dcb, double *rms)
{
	int i,sat;
	char buff[1024],id[32],*label;

	trace(3,"readionexdcb:\n");

	for (i=0;i<MAXSAT;i++) dcb[i]=rms[i]=0.0;

	while (fgets(buff,sizeof(buff),fp)) {
		if (strlen(buff)<60) continue;
		label=buff+60;

		if (strstr(label,"PRN / BIAS / RMS")==label) {

			strncpy(id,buff+3,3); id[3]='\0';

			if (!(sat=satid2no(id))) {
				trace(2,"ionex invalid satellite: %s\n",id);
				continue;
			}
			dcb[sat-1]=str2num(buff, 6,10);
			rms[sat-1]=str2num(buff,16,10);
		}
		else if (strstr(label,"END OF AUX DATA")==label) break;
	}
}

/* read ionex header ---------------------------------------------------------*/
static double readionexh(FILE *fp, double *lats, double *lons, double *hgts,
						 double *rb, double *nexp, double *dcb, double *rms)
{
	double ver=0.0;
	char buff[1024],*label;

	trace(3,"readionexh:\n");

	while (fgets(buff,sizeof(buff),fp)) {

		if (strlen(buff)<60) continue;
		label=buff+60;

		if (strstr(label,"IONEX VERSION / TYPE")==label) {
			if (buff[20]=='I') ver=str2num(buff,0,8);
		}
		else if (strstr(label,"BASE RADIUS")==label) {
			*rb=str2num(buff,0,8);
		}
		else if (strstr(label,"HGT1 / HGT2 / DHGT")==label) {
			hgts[0]=str2num(buff, 2,6);
			hgts[1]=str2num(buff, 8,6);
			hgts[2]=str2num(buff,14,6);
		}
		else if (strstr(label,"LAT1 / LAT2 / DLAT")==label) {
			lats[0]=str2num(buff, 2,6);
			lats[1]=str2num(buff, 8,6);
			lats[2]=str2num(buff,14,6);
		}
		else if (strstr(label,"LON1 / LON2 / DLON")==label) {
			lons[0]=str2num(buff, 2,6);
			lons[1]=str2num(buff, 8,6);
			lons[2]=str2num(buff,14,6);
		}
		else if (strstr(label,"EXPONENT")==label) {
			*nexp=str2num(buff,0,6);
		}
		else if (strstr(label,"START OF AUX DATA")==label&&
			strstr(buff,"DIFFERENTIAL CODE BIASES")) {
				readionexdcb(fp,dcb,rms);
		}
		else if (strstr(label,"END OF HEADER")==label) {
			return ver;
		}
	}
	return 0.0;
}


/* read ionex body -----------------------------------------------------------*/
static int readionexb(FILE *fp, const double *lats, const double *lons,
					  const double *hgts, double rb, double nexp, nav_t *nav)
{
	tec_t *p=NULL;
	gtime_t time={0};
	double lat,lon[3],hgt,x;
	int i,j,k,n,m,index,type=0;
	char buff[1024],*label=buff+60;

	trace(3,"readionexb:\n");

	while (fgets(buff,sizeof(buff),fp)) {

		if (strlen(buff)<60) continue;

		if (strstr(label,"START OF TEC MAP")==label) {
			if ((p=addtec(lats,lons,hgts,rb,nav))) type=1;
		}
		else if (strstr(label,"END OF TEC MAP")==label) {
			type=0;
			p=NULL;
		}
		else if (strstr(label,"START OF RMS MAP")==label) {
			type=2;
			p=NULL;
		}
		else if (strstr(label,"END OF RMS MAP")==label) {
			type=0;
			p=NULL;
		}
		else if (strstr(label,"EPOCH OF CURRENT MAP")==label) {
			if (str2time(buff,0,36,&time)) {
				trace(2,"ionex epoch invalid: %-36.36s\n",buff);
				continue;
			}
			if (type==2) {
				for (i=nav->nt-1;i>=0;i--) {
					if (fabs(timediff(time,nav->tec[i].time))>=1.0) continue;
					p=nav->tec+i;
					break;
				}
			}
			else if (p) p->time=time;
		}
		else if (strstr(label,"LAT/LON1/LON2/DLON/H")==label&&p) {
			lat   =str2num(buff, 2,6);
			lon[0]=str2num(buff, 8,6);
			lon[1]=str2num(buff,14,6);
			lon[2]=str2num(buff,20,6);
			hgt   =str2num(buff,26,6);

			i=getindex(lat,p->lats);
			k=getindex(hgt,p->hgts);
			n=nitem(lon);

			for (m=0;m<n;m++) {
				if (m%16==0&&!fgets(buff,sizeof(buff),fp)) break;

				j=getindex(lon[0]+lon[2]*m,p->lons);
				if ((index=dataindex(i,j,k,p->ndata))<0) continue;

				if ((x=str2num(buff,m%16*5,5))==9999.0) continue;

				if (type==1) p->data[index]=x*pow(10.0,nexp);
				else p->rms[index]=(float)(x*pow(10.0,nexp));
			}
		}
	}
	return 1;
}

/* combine tec grid data -----------------------------------------------------*/
static void combtec(nav_t *nav)
{
	tec_t tmp;
	int i,j,n=0;

	trace(3,"combtec : nav->nt=%d\n",nav->nt);

	for (i=0;i<nav->nt-1;i++) 
	{
		for (j=i+1;j<nav->nt;j++) 
		{
			if (timediff(nav->tec[j].time,nav->tec[i].time)<0.0) 
			{
				tmp=nav->tec[i];
				nav->tec[i]=nav->tec[j];
				nav->tec[j]=tmp;
			}
		}
	}
	for (i=0;i<nav->nt;i++) 
	{
		if (i>0&&timediff(nav->tec[i].time,nav->tec[n-1].time)==0.0) 
		{
			free(nav->tec[n-1].data);
			free(nav->tec[n-1].rms );
			nav->tec[n-1]=nav->tec[i];
			continue;
		}
		nav->tec[n++]=nav->tec[i];
	}
	nav->nt=n;

	trace(4,"combtec : nav->nt=%d\n",nav->nt);
}

/* read ionex tec grid file ----------------------------------------------------
* read ionex ionospheric tec grid file
* args   : char   *file       I   ionex tec grid file
*                                 (wind-card * is expanded)
*          nav_t  *nav        IO  navigation data
*                                 nav->nt, nav->ntmax and nav->tec are modified
*          int    opt         I   read option (1: no clear of tec data,0:clear)
* return : none
* notes  : see ref [1]
*-----------------------------------------------------------------------------*/
extern void readtec(const char *file, nav_t *nav, int opt)
{
	FILE *fp;
	double lats[3]={0},lons[3]={0},hgts[3]={0},rb=0.0,nexp=-1.0;
	double dcb[MAXSAT]={0},rms[MAXSAT]={0};
	int i,n;
	char *efiles[MAXEXFILE];

	trace(3,"readtec : file=%s\n",file);

	/* clear of tec grid data option */
	if (!opt) {
		free(nav->tec); nav->tec=NULL; nav->nt=nav->ntmax=0;
	}
	for (i=0;i<MAXEXFILE;i++) {
		if (!(efiles[i]=(char *)malloc(1024))) {
			for (i--;i>=0;i--) free(efiles[i]);
			return;
		}
	}
	/* expand wild card in file path */
	n=expath(file,efiles,MAXEXFILE);

	for (i=0;i<n;i++) 
	{
		if (!(fp=fopen(efiles[i],"r"))) 
		{
			trace(2,"ionex file open error %s\n",efiles[i]);
			continue;
		}
		/* read ionex header */
		if (readionexh(fp,lats,lons,hgts,&rb,&nexp,dcb,rms)<=0.0) 
		{
			trace(2,"ionex file format error %s\n",efiles[i]);
			continue;
		}
		/* read ionex body */
		readionexb(fp,lats,lons,hgts,rb,nexp,nav);

		fclose(fp);
	}
	for (i=0;i<MAXEXFILE;i++) free(efiles[i]);

	/* combine tec grid data */
	if (nav->nt>0) combtec(nav);

	/* P1-P2 dcb */
	for (i=0;i<MAXSAT;i++) 
	{
		nav->cbias[i][0]=CLIGHT*dcb[i]*1E-9; /* ns->m */
	}
}


/* ionospheric correction ------------------------------------------------------
* compute ionospheric correction
* args   : gtime_t time     I   time
*          nav_t  *nav      I   navigation data
*          int    sat       I   satellite number
*          double *pos      I   receiver position {lat,lon,h} (rad|m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
*          int    ionoopt   I   ionospheric correction option (IONOOPT_???)
*          double *ion      O   ionospheric delay (L1) (m)
*          double *var      O   ionospheric delay (L1) variance (m^2)
* return : status(1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern int ionocorr(gtime_t time, const nav_t *nav, int sat, const double *pos,
					const double *azel, int ionoopt, double *ion, double *var)
{
	trace(4,"ionocorr: time=%s opt=%d sat=%2d pos=%.3f %.3f azel=%.3f %.3f\n",
		time_str(time,3),ionoopt,sat,pos[0]*R2D,pos[1]*R2D,azel[0]*R2D,
		azel[1]*R2D);

	/* broadcast model */
	if (ionoopt==IONOOPT_BRDC) {
		*ion=ionmodel(time,nav->ion_gps,pos,azel);
		*var=SQR(*ion*ERR_BRDCI);
		return 1;
	}

	/* sbas ionosphere model */
	if (ionoopt==IONOOPT_SBAS) 
	{
		return sbsioncorr(time,nav,pos,azel,ion,var);
	}

	/* ionex tec model */
	if (ionoopt==IONOOPT_TEC) 
	{
		return iontec(time,nav,pos,azel,1,ion,var);
	}

	/* qzss broadcast model */
	if (ionoopt==IONOOPT_QZS&&norm(nav->ion_qzs,8)>0.0)
	{
		*ion=ionmodel(time,nav->ion_qzs,pos,azel);
		*var=SQR(*ion*ERR_BRDCI);
		return 1;
	}

	/* lex ionosphere model */
	if (ionoopt==IONOOPT_LEX) 
	{
		return lexioncorr(time,nav,pos,azel,ion,var);
	}
	*ion=0.0;
	*var=ionoopt==IONOOPT_OFF?SQR(ERR_ION):0.0;
	return 1;
}

/* tropospheric correction -----------------------------------------------------
* compute tropospheric correction
* args   : gtime_t time     I   time
*          nav_t  *nav      I   navigation data
*          double *pos      I   receiver position {lat,lon,h} (rad|m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
*          int    tropopt   I   tropospheric correction option (TROPOPT_???)
*          double *trp      O   tropospheric delay (m)
*          double *var      O   tropospheric delay variance (m^2)
* return : status(1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern int tropcorr(gtime_t time, const nav_t *nav, const double *pos,
					const double *azel, int tropopt, double *trp, double *var)
{
	trace(4,"tropcorr: time=%s opt=%d pos=%.3f %.3f azel=%.3f %.3f\n",
		time_str(time,3),tropopt,pos[0]*R2D,pos[1]*R2D,azel[0]*R2D,
		azel[1]*R2D);

	/* saastamoinen model */
	if (tropopt==TROPOPT_SAAS||tropopt==TROPOPT_EST||tropopt==TROPOPT_ESTG) {
		*trp=tropmodel(time,pos,azel,REL_HUMI);
		*var=SQR(ERR_SAAS/(sin(azel[1])+0.1));
		return 1;
	}
	/* sbas troposphere model */
	if (tropopt==TROPOPT_SBAS)
	{
		*trp=sbstropcorr(time,pos,azel,var);
		return 1;
	}
	/* no correction */
	*trp=0.0;
	*var=tropopt==TROPOPT_OFF?SQR(ERR_TROP):0.0;
	return 1;
}