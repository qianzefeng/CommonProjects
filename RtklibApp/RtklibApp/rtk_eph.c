#include "rtk_eph.h"
#include "rtk_trace.h"
#include "rtk_time.h"
#include "rtk_sat.h"
#include "rtk_mat.h"
#include "rtk_sbas.h"
#include "rtk_qzss.h"
#include "rtk_platform.h"
#include "rtk_antenna.h"
#include "rtk_earth.h"


/*------------------------------------------------------------------------------
* ephemeris.c : satellite ephemeris and clock functions
*
*          Copyright (C) 2010-2014 by T.TAKASU, All rights reserved.
*
* references :
*     [1] IS-GPS-200D, Navstar GPS Space Segment/Navigation User Interfaces,
*         7 March, 2006
*     [2] Global Navigation Satellite System GLONASS, Interface Control Document
*         Navigational radiosignal In bands L1, L2, (Edition 5.1), 2008
*     [3] RTCA/DO-229C, Minimum operational performanc standards for global
*         positioning system/wide area augmentation system airborne equipment,
*         RTCA inc, November 28, 2001
*     [4] RTCM Paper, April 12, 2010, Proposed SSR Messages for SV Orbit Clock,
*         Code Biases, URA
*     [5] RTCM Paper 012-2009-SC104-528, January 28, 2009 (previous ver of [4])
*     [6] RTCM Paper 012-2009-SC104-582, February 2, 2010 (previous ver of [4])
*     [7] European GNSS (Galileo) Open Service Signal In Space Interface Control
*         Document, Issue 1, February, 2010
*     [8] Quasi-Zenith Satellite System Navigation Service Interface Control
*         Specification for QZSS (IS-QZSS) V1.1, Japan Aerospace Exploration
*         Agency, July 31, 2009
*     [9] BeiDou navigation satellite system signal in space interface control
*         document open service signal B1I (version 1.0), China Satellite
*         Navigation office, December 2012
*     [10] RTCM Standard 10403.1 - Amendment 5, Differential GNSS (Global
*         Navigation Satellite Systems) Services - version 3, July 1, 2011
*
* version : $Revision:$ $Date:$
* history : 2010/07/28 1.1  moved from rtkcmn.c
*                           added api:
*                               eph2clk(),geph2clk(),seph2clk(),satantoff()
*                               satposs()
*                           changed api:
*                               eph2pos(),geph2pos(),satpos()
*                           deleted api:
*                               satposv(),satposiode()
*           2010/08/26 1.2  add ephemeris option EPHOPT_LEX
*           2010/09/09 1.3  fix problem when precise clock outage
*           2011/01/12 1.4  add api alm2pos()
*                           change api satpos(),satposs()
*                           enable valid unhealthy satellites and output status
*                           fix bug on exception by glonass ephem computation
*           2013/01/10 1.5  support beidou (compass)
*                           use newton's method to solve kepler eq.
*                           update ssr correction algorithm
*           2013/03/20 1.6  fix problem on ssr clock relativitic correction
*           2013/09/01 1.7  support negative pseudorange
*                           fix bug on variance in case of ura ssr = 63
*           2013/11/11 1.8  change constant MAXAGESSR 70.0 -> 90.0
*           2014/10/24 1.9  fix bug on return of var_uraeph() if ura<0||15<ura
*           2014/12/07 1.10 modify MAXDTOE for qzss,gal and bds
*                           test max number of iteration for Kepler
*-----------------------------------------------------------------------------*/

/* constants and macros ------------------------------------------------------*/

#define SQR(x)   ((x)*(x))

#define RE_GLO   6378136.0        /* radius of earth (m)            ref [2] */
#define MU_GPS   3.9860050E14     /* gravitational constant         ref [1] */
#define MU_GLO   3.9860044E14     /* gravitational constant         ref [2] */
#define MU_GAL   3.986004418E14   /* earth gravitational constant   ref [7] */
#define MU_CMP   3.986004418E14   /* earth gravitational constant   ref [9] */
#define J2_GLO   1.0826257E-3     /* 2nd zonal harmonic of geopot   ref [2] */

#define OMGE_GLO 7.292115E-5      /* earth angular velocity (rad/s) ref [2] */
#define OMGE_GAL 7.2921151467E-5  /* earth angular velocity (rad/s) ref [7] */
#define OMGE_CMP 7.292115E-5      /* earth angular velocity (rad/s) ref [9] */

#define SIN_5 -0.0871557427476582 /* sin(-5.0 deg) */
#define COS_5  0.9961946980917456 /* cos(-5.0 deg) */

#define ERREPH_GLO 5.0            /* error of glonass ephemeris (m) */
#define TSTEP    60.0             /* integration step glonass ephemeris (s) */
#define RTOL_KEPLER 1E-14         /* relative tolerance for Kepler equation */

#define DEFURASSR 0.15            /* default accurary of ssr corr (m) */
#define MAXECORSSR 10.0           /* max orbit correction of ssr (m) */
#define MAXCCORSSR (1E-6*CLIGHT)  /* max clock correction of ssr (m) */
#define MAXAGESSR 90.0            /* max age of ssr orbit and clock (s) */
#define MAXAGESSR_HRCLK 10.0      /* max age of ssr high-rate clock (s) */
#define STD_BRDCCLK 30.0          /* error of broadcast clock (m) */

#define MAX_ITER_KEPLER 30        /* max number of iteration of Kelpler */

/* variance by ura ephemeris (ref [1] 20.3.3.3.1.1) --------------------------*/
static double var_uraeph(int ura)
{
	const double ura_value[]=
	{   
		2.4,3.4,4.85,6.85,9.65,13.65,24.0,48.0,96.0,192.0,384.0,768.0,1536.0,
		3072.0,6144.0
	};
	return ura<0||15<ura?SQR(6144.0):SQR(ura_value[ura]);
}

/* variance by ura ssr (ref [4]) ---------------------------------------------*/
static double var_urassr(int ura)
{
	double std;
	if (ura<= 0) return SQR(DEFURASSR);
	if (ura>=63) return SQR(5.4665);
	std=(pow(3.0,(ura>>3)&7)*(1.0+(ura&7)/4.0)-1.0)*1E-3;
	return SQR(std);
}

/* almanac to satellite position and clock bias --------------------------------
* compute satellite position and clock bias with almanac (gps, galileo, qzss)
* args   : gtime_t time     I   time (gpst)
*          alm_t *alm       I   almanac
*          double *rs       O   satellite position (ecef) {x,y,z} (m)
*          double *dts      O   satellite clock bias (s)
* return : none
* notes  : see ref [1],[7],[8]
*-----------------------------------------------------------------------------*/
extern void alm2pos(gtime_t time, const alm_t *alm, double *rs, double *dts)
{
	double tk,M,E,Ek,sinE,cosE,u,r,i,O,x,y,sinO,cosO,cosi,mu;

	trace(4,"alm2pos : time=%s sat=%2d\n",time_str(time,3),alm->sat);

	tk=timediff(time,alm->toa);

	if (alm->A<=0.0) 
	{
		rs[0]=rs[1]=rs[2]=*dts=0.0;
		return;
	}
	mu=satsys(alm->sat,NULL)==SYS_GAL?MU_GAL:MU_GPS;

	M=alm->M0+sqrt(mu/(alm->A*alm->A*alm->A))*tk;
	for (E=M,sinE=Ek=0.0;fabs(E-Ek)>1E-12;)
	{
		Ek=E; sinE=sin(Ek); E=M+alm->e*sinE;
	}
	cosE=cos(E);
	u=atan2(sqrt(1.0-alm->e*alm->e)*sinE,cosE-alm->e)+alm->omg;
	r=alm->A*(1.0-alm->e*cosE);
	i=alm->i0;
	O=alm->OMG0+(alm->OMGd-OMGE)*tk-OMGE*alm->toas;
	x=r*cos(u); y=r*sin(u); sinO=sin(O); cosO=cos(O); cosi=cos(i);
	rs[0]=x*cosO-y*cosi*sinO;
	rs[1]=x*sinO+y*cosi*cosO;
	rs[2]=y*sin(i);
	*dts=alm->f0+alm->f1*tk;
}

/* broadcast ephemeris to satellite clock bias ---------------------------------
* compute satellite clock bias with broadcast ephemeris (gps, galileo, qzss)
* args   : gtime_t time     I   time by satellite clock (gpst)
*          eph_t *eph       I   broadcast ephemeris
* return : satellite clock bias (s) without relativeity correction
* notes  : see ref [1],[7],[8]
*          satellite clock does not include relativity correction and tdg
*-----------------------------------------------------------------------------*/
extern double eph2clk(gtime_t time, const eph_t *eph)
{
	double t;
	int i;

	trace(4,"eph2clk : time=%s sat=%2d\n",time_str(time,3),eph->sat);

	t=timediff(time,eph->toc);

	for (i=0;i<2;i++) 
	{
		t-=eph->f0+eph->f1*t+eph->f2*t*t;
	}
	return eph->f0+eph->f1*t+eph->f2*t*t;
}

/* broadcast ephemeris to satellite position and clock bias --------------------
* compute satellite position and clock bias with broadcast ephemeris (gps,
* galileo, qzss)
* args   : gtime_t time     I   time (gpst)
*          eph_t *eph       I   broadcast ephemeris
*          double *rs       O   satellite position (ecef) {x,y,z} (m)
*          double *dts      O   satellite clock bias (s)
*          double *var      O   satellite position and clock variance (m^2)
* return : none
* notes  : see ref [1],[7],[8]
*          satellite clock includes relativity correction without code bias
*          (tgd or bgd)
*-----------------------------------------------------------------------------*/
extern void eph2pos(gtime_t time, const eph_t *eph, double *rs, double *dts,
					double *var)
{
	double tk,M,E,Ek,sinE,cosE,u,r,i,O,sin2u,cos2u,x,y,sinO,cosO,cosi,mu,omge;
	double xg,yg,zg,sino,coso;
	int n,sys,prn;

	trace(4,"eph2pos : time=%s sat=%2d\n",time_str(time,3),eph->sat);

	if (eph->A<=0.0) 
	{
		rs[0]=rs[1]=rs[2]=*dts=*var=0.0;
		return;
	}
	tk=timediff(time,eph->toe);

	switch ((sys=satsys(eph->sat,&prn))) 
	{
		case SYS_GAL: mu=MU_GAL; omge=OMGE_GAL; break;
		case SYS_CMP: mu=MU_CMP; omge=OMGE_CMP; break;
		default:      mu=MU_GPS; omge=OMGE;     break;
	}
	M=eph->M0+(sqrt(mu/(eph->A*eph->A*eph->A))+eph->deln)*tk;

	for (n=0,E=M,Ek=0.0;fabs(E-Ek)>RTOL_KEPLER&&n<MAX_ITER_KEPLER;n++) 
	{
		Ek=E; E-=(E-eph->e*sin(E)-M)/(1.0-eph->e*cos(E));
	}
	if (n>=MAX_ITER_KEPLER) 
	{
		trace(2,"kepler iteration overflow sat=%2d\n",eph->sat);
		return;
	}
	sinE=sin(E); cosE=cos(E);

	trace(4,"kepler: sat=%2d e=%8.5f n=%2d del=%10.3e\n",eph->sat,eph->e,n,E-Ek);

	u=atan2(sqrt(1.0-eph->e*eph->e)*sinE,cosE-eph->e)+eph->omg;
	r=eph->A*(1.0-eph->e*cosE);
	i=eph->i0+eph->idot*tk;
	sin2u=sin(2.0*u); cos2u=cos(2.0*u);
	u+=eph->cus*sin2u+eph->cuc*cos2u;
	r+=eph->crs*sin2u+eph->crc*cos2u;
	i+=eph->cis*sin2u+eph->cic*cos2u;
	x=r*cos(u); y=r*sin(u); cosi=cos(i);

	/* beidou geo satellite (ref [9]) */
	if (sys==SYS_CMP&&prn<=5) 
	{
		O=eph->OMG0+eph->OMGd*tk-omge*eph->toes;
		sinO=sin(O); cosO=cos(O);
		xg=x*cosO-y*cosi*sinO;
		yg=x*sinO+y*cosi*cosO;
		zg=y*sin(i);
		sino=sin(omge*tk); coso=cos(omge*tk);
		rs[0]= xg*coso+yg*sino*COS_5+zg*sino*SIN_5;
		rs[1]=-xg*sino+yg*coso*COS_5+zg*coso*SIN_5;
		rs[2]=-yg*SIN_5+zg*COS_5;
	}
	else
	{
		O=eph->OMG0+(eph->OMGd-omge)*tk-omge*eph->toes;
		sinO=sin(O); cosO=cos(O);
		rs[0]=x*cosO-y*cosi*sinO;
		rs[1]=x*sinO+y*cosi*cosO;
		rs[2]=y*sin(i);
	}
	tk=timediff(time,eph->toc);
	*dts=eph->f0+eph->f1*tk+eph->f2*tk*tk;

	/* relativity correction */
	*dts-=2.0*sqrt(mu*eph->A)*eph->e*sinE/SQR(CLIGHT);

	/* position and clock error variance */
	*var=var_uraeph(eph->sva);
}

/* glonass orbit differential equations --------------------------------------*/
static void deq(const double *x, double *xdot, const double *acc)
{
	double a,b,c,r2=dot(x,x,3),r3=r2*sqrt(r2),omg2=SQR(OMGE_GLO);

	if (r2<=0.0) 
	{
		xdot[0]=xdot[1]=xdot[2]=xdot[3]=xdot[4]=xdot[5]=0.0;
		return;
	}
	/* ref [2] A.3.1.2 with bug fix for xdot[4],xdot[5] */
	a=1.5*J2_GLO*MU_GLO*SQR(RE_GLO)/r2/r3; /* 3/2*J2*mu*Ae^2/r^5 */
	b=5.0*x[2]*x[2]/r2;                    /* 5*z^2/r^2 */
	c=-MU_GLO/r3-a*(1.0-b);                /* -mu/r^3-a(1-b) */
	xdot[0]=x[3]; xdot[1]=x[4]; xdot[2]=x[5];
	xdot[3]=(c+omg2)*x[0]+2.0*OMGE_GLO*x[4]+acc[0];
	xdot[4]=(c+omg2)*x[1]-2.0*OMGE_GLO*x[3]+acc[1];
	xdot[5]=(c-2.0*a)*x[2]+acc[2];
}

/* glonass position and velocity by numerical integration --------------------*/
static void glorbit(double t, double *x, const double *acc)
{
	double k1[6],k2[6],k3[6],k4[6],w[6];
	int i;

	deq(x,k1,acc); for (i=0;i<6;i++) w[i]=x[i]+k1[i]*t/2.0;
	deq(w,k2,acc); for (i=0;i<6;i++) w[i]=x[i]+k2[i]*t/2.0;
	deq(w,k3,acc); for (i=0;i<6;i++) w[i]=x[i]+k3[i]*t;
	deq(w,k4,acc);
	for (i=0;i<6;i++) x[i]+=(k1[i]+2.0*k2[i]+2.0*k3[i]+k4[i])*t/6.0;
}

/* glonass ephemeris to satellite clock bias -----------------------------------
* compute satellite clock bias with glonass ephemeris
* args   : gtime_t time     I   time by satellite clock (gpst)
*          geph_t *geph     I   glonass ephemeris
* return : satellite clock bias (s)
* notes  : see ref [2]
*-----------------------------------------------------------------------------*/
extern double geph2clk(gtime_t time, const geph_t *geph)
{
	double t;
	int i;

	trace(4,"geph2clk: time=%s sat=%2d\n",time_str(time,3),geph->sat);

	t=timediff(time,geph->toe);

	for (i=0;i<2;i++) 
	{
		t-=-geph->taun+geph->gamn*t;
	}
	return -geph->taun+geph->gamn*t;
}

/* glonass ephemeris to satellite position and clock bias ----------------------
* compute satellite position and clock bias with glonass ephemeris
* args   : gtime_t time     I   time (gpst)
*          geph_t *geph     I   glonass ephemeris
*          double *rs       O   satellite position {x,y,z} (ecef) (m)
*          double *dts      O   satellite clock bias (s)
*          double *var      O   satellite position and clock variance (m^2)
* return : none
* notes  : see ref [2]
*-----------------------------------------------------------------------------*/
extern void geph2pos(gtime_t time, const geph_t *geph, double *rs, double *dts,
					 double *var)
{
	double t,tt,x[6];
	int i;

	trace(4,"geph2pos: time=%s sat=%2d\n",time_str(time,3),geph->sat);

	t=timediff(time,geph->toe);

	*dts=-geph->taun+geph->gamn*t;

	for (i=0;i<3;i++) 
	{
		x[i  ]=geph->pos[i];
		x[i+3]=geph->vel[i];
	}
	for (tt=t<0.0?-TSTEP:TSTEP;fabs(t)>1E-9;t-=tt)
	{
		if (fabs(t)<TSTEP) tt=t;
		glorbit(tt,x,geph->acc);
	}
	for (i=0;i<3;i++) rs[i]=x[i];

	*var=SQR(ERREPH_GLO);
}

/* sbas ephemeris to satellite clock bias --------------------------------------
* compute satellite clock bias with sbas ephemeris
* args   : gtime_t time     I   time by satellite clock (gpst)
*          seph_t *seph     I   sbas ephemeris
* return : satellite clock bias (s)
* notes  : see ref [3]
*-----------------------------------------------------------------------------*/
extern double seph2clk(gtime_t time, const seph_t *seph)
{
	double t;
	int i;

	trace(4,"seph2clk: time=%s sat=%2d\n",time_str(time,3),seph->sat);

	t=timediff(time,seph->t0);

	for (i=0;i<2;i++) 
	{
		t-=seph->af0+seph->af1*t;
	}
	return seph->af0+seph->af1*t;
}

/* sbas ephemeris to satellite position and clock bias -------------------------
* compute satellite position and clock bias with sbas ephemeris
* args   : gtime_t time     I   time (gpst)
*          seph_t  *seph    I   sbas ephemeris
*          double  *rs      O   satellite position {x,y,z} (ecef) (m)
*          double  *dts     O   satellite clock bias (s)
*          double  *var     O   satellite position and clock variance (m^2)
* return : none
* notes  : see ref [3]
*-----------------------------------------------------------------------------*/
extern void seph2pos(gtime_t time, const seph_t *seph, double *rs, double *dts,
					 double *var)
{
	double t;
	int i;

	trace(4,"seph2pos: time=%s sat=%2d\n",time_str(time,3),seph->sat);

	t=timediff(time,seph->t0);

	for (i=0;i<3;i++) 
	{
		rs[i]=seph->pos[i]+seph->vel[i]*t+seph->acc[i]*t*t/2.0;
	}
	*dts=seph->af0+seph->af1*t;

	*var=var_uraeph(seph->sva);
}

/* select ephememeris --------------------------------------------------------*/
static eph_t *seleph(gtime_t time, int sat, int iode, const nav_t *nav)
{
	double t,tmax,tmin;
	int i,j=-1;

	trace(4,"seleph  : time=%s sat=%2d iode=%d\n",time_str(time,3),sat,iode);

	switch (satsys(sat,NULL))
	{
		case SYS_QZS: tmax=MAXDTOE_QZS+1.0; break;
		case SYS_GAL: tmax=MAXDTOE_GAL+1.0; break;
		case SYS_CMP: tmax=MAXDTOE_CMP+1.0; break;
		default: tmax=MAXDTOE+1.0; break;
	}
	tmin=tmax+1.0;

	for (i=0;i<nav->n;i++) 
	{
		if (nav->eph[i].sat!=sat) continue;
		if (iode>=0&&nav->eph[i].iode!=iode) continue;
		if ((t=fabs(timediff(nav->eph[i].toe,time)))>tmax) continue;
		if (iode>=0) return nav->eph+i;
		if (t<=tmin) {j=i; tmin=t;} /* toe closest to time */
	}
	if (iode>=0||j<0) 
	{
		trace(2,"no broadcast ephemeris: %s sat=%2d iode=%3d\n",time_str(time,0),
			sat,iode);
		return NULL;
	}
	return nav->eph+j;
}

/* select glonass ephememeris ------------------------------------------------*/
static geph_t *selgeph(gtime_t time, int sat, int iode, const nav_t *nav)
{
	double t,tmax=MAXDTOE_GLO,tmin=tmax+1.0;
	int i,j=-1;

	trace(4,"selgeph : time=%s sat=%2d iode=%2d\n",time_str(time,3),sat,iode);

	for (i=0;i<nav->ng;i++) 
	{
		if (nav->geph[i].sat!=sat) continue;
		if (iode>=0&&nav->geph[i].iode!=iode) continue;
		if ((t=fabs(timediff(nav->geph[i].toe,time)))>tmax) continue;
		if (iode>=0) return nav->geph+i;
		if (t<=tmin) {j=i; tmin=t;} /* toe closest to time */
	}
	if (iode>=0||j<0) 
	{
		trace(3,"no glonass ephemeris  : %s sat=%2d iode=%2d\n",time_str(time,0),
			sat,iode);
		return NULL;
	}
	return nav->geph+j;
}

/* select sbas ephememeris ---------------------------------------------------*/
static seph_t *selseph(gtime_t time, int sat, const nav_t *nav)
{
	double t,tmax=MAXDTOE_SBS,tmin=tmax+1.0;
	int i,j=-1;

	trace(4,"selseph : time=%s sat=%2d\n",time_str(time,3),sat);

	for (i=0;i<nav->ns;i++) 
	{
		if (nav->seph[i].sat!=sat) continue;
		if ((t=fabs(timediff(nav->seph[i].t0,time)))>tmax) continue;
		if (t<=tmin) {j=i; tmin=t;} /* toe closest to time */
	}
	if (j<0) 
	{
		trace(3,"no sbas ephemeris     : %s sat=%2d\n",time_str(time,0),sat);
		return NULL;
	}
	return nav->seph+j;
}

/* satellite clock with broadcast ephemeris ----------------------------------*/
static int ephclk(gtime_t time, gtime_t teph, int sat, const nav_t *nav,
				  double *dts)
{
	eph_t  *eph;
	geph_t *geph;
	seph_t *seph;
	int sys;

	trace(4,"ephclk  : time=%s sat=%2d\n",time_str(time,3),sat);

	sys=satsys(sat,NULL);

	if (sys==SYS_GPS||sys==SYS_GAL||sys==SYS_QZS||sys==SYS_CMP) 
	{
		if (!(eph=seleph(teph,sat,-1,nav))) return 0;
		*dts=eph2clk(time,eph);
	}
	else if (sys==SYS_GLO)
	{
		if (!(geph=selgeph(teph,sat,-1,nav))) return 0;
		*dts=geph2clk(time,geph);
	}
	else if (sys==SYS_SBS)
	{
		if (!(seph=selseph(teph,sat,nav))) return 0;
		*dts=seph2clk(time,seph);
	}
	else return 0;

	return 1;
}

/* satellite position and clock by broadcast ephemeris -----------------------*/
static int ephpos(gtime_t time, gtime_t teph, int sat, const nav_t *nav,
				  int iode, double *rs, double *dts, double *var, int *svh)
{
	eph_t  *eph;
	geph_t *geph;
	seph_t *seph;
	double rst[3],dtst[1],tt=1E-3;
	int i,sys;

	trace(4,"ephpos  : time=%s sat=%2d iode=%d\n",time_str(time,3),sat,iode);

	sys=satsys(sat,NULL);

	*svh=-1;

	if (sys==SYS_GPS||sys==SYS_GAL||sys==SYS_QZS||sys==SYS_CMP) 
	{
		if (!(eph=seleph(teph,sat,iode,nav))) return 0;

		eph2pos(time,eph,rs,dts,var);
		time=timeadd(time,tt);
		eph2pos(time,eph,rst,dtst,var);
		*svh=eph->svh;
	}
	else if (sys==SYS_GLO) 
	{
		if (!(geph=selgeph(teph,sat,iode,nav))) return 0;
		geph2pos(time,geph,rs,dts,var);
		time=timeadd(time,tt);
		geph2pos(time,geph,rst,dtst,var);
		*svh=geph->svh;
	}
	else if (sys==SYS_SBS) 
	{
		if (!(seph=selseph(teph,sat,nav))) return 0;

		seph2pos(time,seph,rs,dts,var);
		time=timeadd(time,tt);
		seph2pos(time,seph,rst,dtst,var);
		*svh=seph->svh;
	}
	else return 0;

	/* satellite velocity and clock drift by differential approx */
	for (i=0;i<3;i++) rs[i+3]=(rst[i]-rs[i])/tt;
	dts[1]=(dtst[0]-dts[0])/tt;

	return 1;
}

/* satellite position and clock with sbas correction -------------------------*/
static int satpos_sbas(gtime_t time, gtime_t teph, int sat, const nav_t *nav,
					   double *rs, double *dts, double *var, int *svh)
{
	const sbssatp_t *sbs;
	int i;

	trace(4,"satpos_sbas: time=%s sat=%2d\n",time_str(time,3),sat);

	/* search sbas satellite correciton */
	for (i=0;i<nav->sbssat.nsat;i++)
	{
		sbs=nav->sbssat.sat+i;
		if (sbs->sat==sat) break;
	}

	if (i>=nav->sbssat.nsat) 
	{
		trace(2,"no sbas correction for orbit: %s sat=%2d\n",time_str(time,0),sat);
		ephpos(time,teph,sat,nav,-1,rs,dts,var,svh);
		*svh=-1;
		return 0;
	}
	/* satellite postion and clock by broadcast ephemeris */
	if (!ephpos(time,teph,sat,nav,sbs->lcorr.iode,rs,dts,var,svh)) return 0;

	/* sbas satellite correction (long term and fast) */
	if (sbssatcorr(time,sat,nav,rs,dts,var)) return 1;
	*svh=-1;
	return 0;
}

/* satellite position and clock with ssr correction --------------------------*/
static int satpos_ssr(gtime_t time, gtime_t teph, int sat, const nav_t *nav,
					  int opt, double *rs, double *dts, double *var, int *svh)
{
	const ssr_t *ssr;
	eph_t *eph;
	double t1,t2,t3,er[3],ea[3],ec[3],rc[3],deph[3],dclk,dant[3]={0},tk;
	int i,sys;

	trace(4,"satpos_ssr: time=%s sat=%2d\n",time_str(time,3),sat);

	ssr=nav->ssr+sat-1;

	if (!ssr->t0[0].time) 
	{
		trace(2,"no ssr orbit correction: %s sat=%2d\n",time_str(time,0),sat);
		return 0;
	}
	if (!ssr->t0[1].time) 
	{
		trace(2,"no ssr clock correction: %s sat=%2d\n",time_str(time,0),sat);
		return 0;
	}
	/* inconsistency between orbit and clock correction */
	if (ssr->iod[0]!=ssr->iod[1]) 
	{
		trace(2,"inconsist ssr correction: %s sat=%2d iod=%d %d\n",
			time_str(time,0),sat,ssr->iod[0],ssr->iod[1]);
		*svh=-1;
		return 0;
	}
	t1=timediff(time,ssr->t0[0]);
	t2=timediff(time,ssr->t0[1]);
	t3=timediff(time,ssr->t0[2]);

	/* ssr orbit and clock correction (ref [4]) */
	if (fabs(t1)>MAXAGESSR||fabs(t2)>MAXAGESSR) 
	{
		trace(2,"age of ssr error: %s sat=%2d t=%.0f %.0f\n",time_str(time,0),
			sat,t1,t2);
		*svh=-1;
		return 0;
	}
	if (ssr->udi[0]>=1.0) t1-=ssr->udi[0]/2.0;
	if (ssr->udi[1]>=1.0) t2-=ssr->udi[0]/2.0;

	for (i=0;i<3;i++) deph[i]=ssr->deph[i]+ssr->ddeph[i]*t1;
	dclk=ssr->dclk[0]+ssr->dclk[1]*t2+ssr->dclk[2]*t2*t2;

	/* ssr highrate clock correction (ref [4]) */
	if (ssr->iod[0]==ssr->iod[2]&&ssr->t0[2].time&&fabs(t3)<MAXAGESSR_HRCLK) 
	{
		dclk+=ssr->hrclk;
	}
	if (norm(deph,3)>MAXECORSSR||fabs(dclk)>MAXCCORSSR) 
	{
		trace(3,"invalid ssr correction: %s deph=%.1f dclk=%.1f\n",
			time_str(time,0),norm(deph,3),dclk);
		*svh=-1;
		return 0;
	}
	/* satellite postion and clock by broadcast ephemeris */
	if (!ephpos(time,teph,sat,nav,ssr->iode,rs,dts,var,svh)) return 0;

	/* satellite clock for gps, galileo and qzss */
	sys=satsys(sat,NULL);
	if (sys==SYS_GPS||sys==SYS_GAL||sys==SYS_QZS||sys==SYS_CMP) 
	{
		if (!(eph=seleph(teph,sat,ssr->iode,nav))) return 0;

		/* satellite clock by clock parameters */
		tk=timediff(time,eph->toc);
		dts[0]=eph->f0+eph->f1*tk+eph->f2*tk*tk;
		dts[1]=eph->f1+2.0*eph->f2*tk;

		/* relativity correction */
		dts[0]-=2.0*dot(rs,rs+3,3)/CLIGHT/CLIGHT;
	}
	/* radial-along-cross directions in ecef */
	if (!normv3(rs+3,ea)) return 0;
	cross3(rs,rs+3,rc);
	if (!normv3(rc,ec)) 
	{
		*svh=-1;
		return 0;
	}
	cross3(ea,ec,er);

	/* satellite antenna offset correction */
	if (opt) 
	{
		satantoff(time,rs,sat,nav,dant);
	}
	for (i=0;i<3;i++) {
		rs[i]+=-(er[i]*deph[0]+ea[i]*deph[1]+ec[i]*deph[2])+dant[i];
	}
	/* t_corr = t_sv - (dts(brdc) + dclk(ssr) / CLIGHT) (ref [10] eq.3.12-7) */
	dts[0]+=dclk/CLIGHT;

	/* variance by ssr ura */
	*var=var_urassr(ssr->ura);

	trace(5,"satpos_ssr: %s sat=%2d deph=%6.3f %6.3f %6.3f er=%6.3f %6.3f %6.3f dclk=%6.3f var=%6.3f\n",
		time_str(time,2),sat,deph[0],deph[1],deph[2],er[0],er[1],er[2],dclk,*var);

	return 1;
}

/* satellite position and clock ------------------------------------------------
* compute satellite position, velocity and clock
* args   : gtime_t time     I   time (gpst)
*          gtime_t teph     I   time to select ephemeris (gpst)
*          int    sat       I   satellite number
*          nav_t  *nav      I   navigation data
*          int    ephopt    I   ephemeris option (EPHOPT_???)
*          double *rs       O   sat position and velocity (ecef)
*                               {x,y,z,vx,vy,vz} (m|m/s)
*          double *dts      O   sat clock {bias,drift} (s|s/s)
*          double *var      O   sat position and clock error variance (m^2)
*          int    *svh      O   sat health flag (-1:correction not available)
* return : status (1:ok,0:error)
* notes  : satellite position is referenced to antenna phase center
*          satellite clock does not include code bias correction (tgd or bgd)
*-----------------------------------------------------------------------------*/
extern int satpos(gtime_t time, gtime_t teph, int sat, int ephopt,
				  const nav_t *nav, double *rs, double *dts, double *var,
				  int *svh)
{
	trace(4,"satpos  : time=%s sat=%2d ephopt=%d\n",time_str(time,3),sat,ephopt);

	*svh=0;

	switch (ephopt) 
	{
		case EPHOPT_BRDC  : return ephpos     (time,teph,sat,nav,-1,rs,dts,var,svh);
		case EPHOPT_SBAS  : return satpos_sbas(time,teph,sat,nav,   rs,dts,var,svh);
		case EPHOPT_SSRAPC: return satpos_ssr (time,teph,sat,nav, 0,rs,dts,var,svh);
		case EPHOPT_SSRCOM: return satpos_ssr (time,teph,sat,nav, 1,rs,dts,var,svh);
		case EPHOPT_PREC  :
			if (!peph2pos(time,sat,nav,1,rs,dts,var)) break; else return 1;
		case EPHOPT_LEX   :
			if (!lexeph2pos(time,sat,nav,rs,dts,var)) break; else return 1;
	}
	*svh=-1;
	return 0;
}

/* satellite positions and clocks ----------------------------------------------
* compute satellite positions, velocities and clocks
* args   : gtime_t teph     I   time to select ephemeris (gpst)
*          obsd_t *obs      I   observation data
*          int    n         I   number of observation data
*          nav_t  *nav      I   navigation data
*          int    ephopt    I   ephemeris option (EPHOPT_???)
*          double *rs       O   satellite positions and velocities (ecef)
*          double *dts      O   satellite clocks
*          double *var      O   sat position and clock error variances (m^2)
*          int    *svh      O   sat health flag (-1:correction not available)
* return : none
* notes  : rs [(0:2)+i*6]= obs[i] sat position {x,y,z} (m)
*          rs [(3:5)+i*6]= obs[i] sat velocity {vx,vy,vz} (m/s)
*          dts[(0:1)+i*2]= obs[i] sat clock {bias,drift} (s|s/s)
*          var[i]        = obs[i] sat position and clock error variance (m^2)
*          svh[i]        = obs[i] sat health flag
*          if no navigation data, set 0 to rs[], dts[], var[] and svh[]
*          satellite position and clock are values at signal transmission time
*          satellite position is referenced to antenna phase center
*          satellite clock does not include code bias correction (tgd or bgd)
*          any pseudorange and broadcast ephemeris are always needed to get
*          signal transmission time
*-----------------------------------------------------------------------------*/
extern void satposs(gtime_t teph, const obsd_t *obs, int n, const nav_t *nav,
					int ephopt, double *rs, double *dts, double *var, int *svh)
{
	gtime_t time[MAXOBS]={{0}};
	double dt,pr;
	int i,j;

	trace(3,"satposs : teph=%s n=%d ephopt=%d\n",time_str(teph,3),n,ephopt);

	for (i=0;i<n&&i<MAXOBS;i++) 
	{
		for (j=0;j<6;j++) rs [j+i*6]=0.0;
		for (j=0;j<2;j++) dts[j+i*2]=0.0;
		var[i]=0.0; svh[i]=0;

		/* search any psuedorange */
		for (j=0,pr=0.0;j<NFREQ;j++) if ((pr=obs[i].P[j])!=0.0) break;

		if (j>=NFREQ) 
		{
			trace(2,"no pseudorange %s sat=%2d\n",time_str(obs[i].time,3),obs[i].sat);
			continue;
		}
		/* transmission time by satellite clock */
		time[i]=timeadd(obs[i].time,-pr/CLIGHT);

		/* satellite clock bias by broadcast ephemeris */
		if (!ephclk(time[i],teph,obs[i].sat,nav,&dt)) 
		{
			trace(2,"no broadcast clock %s sat=%2d\n",time_str(time[i],3),obs[i].sat);
			continue;
		}
		time[i]=timeadd(time[i],-dt);

		/* satellite position and clock at transmission time */
		if (!satpos(time[i],teph,obs[i].sat,ephopt,nav,rs+i*6,dts+i*2,var+i,
			svh+i)) 
		{
			trace(3,"no ephemeris %s sat=%2d\n",time_str(time[i],3),obs[i].sat);
			continue;
		}
		/* if no precise clock available, use broadcast clock instead */
		if (dts[i*2]==0.0) 
		{
			if (!ephclk(time[i],teph,obs[i].sat,nav,dts+i*2)) continue;
			dts[1+i*2]=0.0;
			*var=SQR(STD_BRDCCLK);
		}
	}
	for (i=0;i<n&&i<MAXOBS;i++) 
	{
		trace(4,"%s sat=%2d rs=%13.3f %13.3f %13.3f dts=%12.3f var=%7.3f svh=%02X\n",
			time_str(time[i],6),obs[i].sat,rs[i*6],rs[1+i*6],rs[2+i*6],
			dts[i*2]*1E9,var[i],svh[i]);
	}
}



/*------------------------------------------------------------------------------
* preceph.c : precise ephemeris and clock functions
*
*          Copyright (C) 2007-2013 by T.TAKASU, All rights reserved.
*
* references :
*     [1] S.Hilla, The Extended Standard Product 3 Orbit Format (SP3-c),
*         12 February, 2007
*     [2] J.Ray, W.Gurtner, RINEX Extensions to Handle Clock Information,
*         27 August, 1998
*     [3] D.D.McCarthy, IERS Technical Note 21, IERS Conventions 1996, July 1996
*     [4] D.A.Vallado, Fundamentals of Astrodynamics and Applications 2nd ed,
*         Space Technology Library, 2004
*
* version : $Revision: 1.1 $ $Date: 2008/07/17 21:48:06 $
* history : 2009/01/18 1.0  new
*           2009/01/31 1.1  fix bug on numerical error to read sp3a ephemeris
*           2009/05/15 1.2  support glonass,galileo,qzs
*           2009/12/11 1.3  support wild-card expansion of file path
*           2010/07/21 1.4  added api:
*                               eci2ecef(),sunmoonpos(),peph2pos(),satantoff(),
*                               readdcb()
*                           changed api:
*                               readsp3()
*                           deleted api:
*                               eph2posp()
*           2010/09/09 1.5  fix problem when precise clock outage
*           2011/01/23 1.6  support qzss satellite code
*           2011/09/12 1.7  fix problem on precise clock outage
*                           move sunmmonpos() to rtkcmn.c
*           2011/12/01 1.8  modify api readsp3()
*                           precede later ephemeris if ephemeris is NULL
*                           move eci2ecef() to rtkcmn.c
*           2013/05/08 1.9  fix bug on computing std-dev of precise clocks
*           2013/11/20 1.10 modify option for api readsp3()
*           2014/04/03 1.11 accept extenstion including sp3,eph,SP3,EPH
*           2014/05/23 1.12 add function to read sp3 velocity records
*                           change api: satantoff()
*           2014/08/31 1.13 add member cov and vco in peph_t sturct
*           2014/10/13 1.14 fix bug on clock error variance in peph2pos()
*           2015/05/10 1.15 add api readfcb()
*                           modify api readdcb()
*-----------------------------------------------------------------------------*/


#define NMAX        10              /* order of polynomial interpolation */
#define MAXDTE      900.0           /* max time difference to ephem time (s) */
#define EXTERR_CLK  1E-3            /* extrapolation error for clock (m/s) */
#define EXTERR_EPH  5E-7            /* extrapolation error for ephem (m/s^2) */

/* satellite code to satellite system ----------------------------------------*/
static int code2sys(char code)
{
	if (code=='G'||code==' ') return SYS_GPS;
	if (code=='R') return SYS_GLO;
	if (code=='E') return SYS_GAL; /* extension to sp3-c */
	if (code=='J') return SYS_QZS; /* extension to sp3-c */
	if (code=='C') return SYS_CMP; /* extension to sp3-c */
	if (code=='L') return SYS_LEO; /* extension to sp3-c */
	return SYS_NONE;
}
/* read sp3 header -----------------------------------------------------------*/
static int readsp3h(FILE *fp, gtime_t *time, char *type, int *sats,
					double *bfact, char *tsys)
{
	int i,j,k=0,ns=0,sys,prn;
	char buff[1024];

	trace(3,"readsp3h:\n");

	for (i=0;i<22;i++) {
		if (!fgets(buff,sizeof(buff),fp)) break;

		if (i==0) {
			*type=buff[2];
			if (str2time(buff,3,28,time)) return 0;
		}
		else if (2<=i&&i<=6) {
			if (i==2) {
				ns=(int)str2num(buff,4,2);
			}
			for (j=0;j<17&&k<ns;j++) {
				sys=code2sys(buff[9+3*j]);
				prn=(int)str2num(buff,10+3*j,2);
				if (k<MAXSAT) sats[k++]=satno(sys,prn);
			}
		}
		else if (i==12) {
			strncpy(tsys,buff+9,3); tsys[3]='\0';
		}
		else if (i==14) {
			bfact[0]=str2num(buff, 3,10);
			bfact[1]=str2num(buff,14,12);
		}
	}
	return ns;
}
/* add precise ephemeris -----------------------------------------------------*/
static int addpeph(nav_t *nav, peph_t *peph)
{
	peph_t *nav_peph;

	if (nav->ne>=nav->nemax) {
		nav->nemax+=256;
		if (!(nav_peph=(peph_t *)realloc(nav->peph,sizeof(peph_t)*nav->nemax))) {
			trace(1,"readsp3b malloc error n=%d\n",nav->nemax);
			free(nav->peph); nav->peph=NULL; nav->ne=nav->nemax=0;
			return 0;
		}
		nav->peph=nav_peph;
	}
	nav->peph[nav->ne++]=*peph;
	return 1;
}
/* read sp3 body -------------------------------------------------------------*/
static void readsp3b(FILE *fp, char type, int *sats, int ns, double *bfact,
					 char *tsys, int index, int opt, nav_t *nav)
{
	peph_t peph;
	gtime_t time;
	double val,std,base;
	int i,j,sat,sys,prn,n=ns*(type=='P'?1:2),pred_o,pred_c,v;
	char buff[1024];

	trace(3,"readsp3b: type=%c ns=%d index=%d opt=%d\n",type,ns,index,opt);

	while (fgets(buff,sizeof(buff),fp)) {

		if (!strncmp(buff,"EOF",3)) break;

		if (buff[0]!='*'||str2time(buff,3,28,&time)) {
			trace(2,"sp3 invalid epoch %31.31s\n",buff);
			continue;
		}
		if (!strcmp(tsys,"UTC")) time=utc2gpst(time); /* utc->gpst */
		peph.time =time;
		peph.index=index;

		for (i=0;i<MAXSAT;i++) {
			for (j=0;j<4;j++) {
				peph.pos[i][j]=0.0;
				peph.std[i][j]=0.0f;
				peph.vel[i][j]=0.0;
				peph.vst[i][j]=0.0f;
			}
			for (j=0;j<3;j++) {
				peph.cov[i][j]=0.0f;
				peph.vco[i][j]=0.0f;
			}
		}
		for (i=pred_o=pred_c=v=0;i<n&&fgets(buff,sizeof(buff),fp);i++) {

			if (strlen(buff)<4||(buff[0]!='P'&&buff[0]!='V')) continue;

			sys=buff[1]==' '?SYS_GPS:code2sys(buff[1]);
			prn=(int)str2num(buff,2,2);
			if      (sys==SYS_SBS) prn+=100;
			else if (sys==SYS_QZS) prn+=192; /* extension to sp3-c */

			if (!(sat=satno(sys,prn))) continue;

			if (buff[0]=='P') {
				pred_c=strlen(buff)>=76&&buff[75]=='P';
				pred_o=strlen(buff)>=80&&buff[79]=='P';
			}
			for (j=0;j<4;j++) {

				/* read option for predicted value */
				if (j< 3&&(opt&1)&& pred_o) continue;
				if (j< 3&&(opt&2)&&!pred_o) continue;
				if (j==3&&(opt&1)&& pred_c) continue;
				if (j==3&&(opt&2)&&!pred_c) continue;

				val=str2num(buff, 4+j*14,14);
				std=str2num(buff,61+j* 3,j<3?2:3);

				if (buff[0]=='P') { /* position */
					if (val!=0.0&&fabs(val-999999.999999)>=1E-6) {
						peph.pos[sat-1][j]=val*(j<3?1000.0:1E-6);
						v=1; /* valid epoch */
					}
					if ((base=bfact[j<3?0:1])>0.0&&std>0.0) {
						peph.std[sat-1][j]=(float)(pow(base,std)*(j<3?1E-3:1E-12));
					}
				}
				else if (v) { /* velocity */
					if (val!=0.0&&fabs(val-999999.999999)>=1E-6) {
						peph.vel[sat-1][j]=val*(j<3?0.1:1E-10);
					}
					if ((base=bfact[j<3?0:1])>0.0&&std>0.0) {
						peph.vst[sat-1][j]=(float)(pow(base,std)*(j<3?1E-7:1E-16));
					}
				}
			}
		}
		if (v) {
			if (!addpeph(nav,&peph)) return;
		}
	}
}
/* compare precise ephemeris -------------------------------------------------*/
static int cmppeph(const void *p1, const void *p2)
{
	peph_t *q1=(peph_t *)p1,*q2=(peph_t *)p2;
	double tt=timediff(q1->time,q2->time);
	return tt<-1E-9?-1:(tt>1E-9?1:q1->index-q2->index);
}
/* combine precise ephemeris -------------------------------------------------*/
static void combpeph(nav_t *nav, int opt)
{
	int i,j,k,m;

	trace(3,"combpeph: ne=%d\n",nav->ne);

	qsort(nav->peph,nav->ne,sizeof(peph_t),cmppeph);

	if (opt&4) return;

	for (i=0,j=1;j<nav->ne;j++) {

		if (fabs(timediff(nav->peph[i].time,nav->peph[j].time))<1E-9) {

			for (k=0;k<MAXSAT;k++) {
				if (norm(nav->peph[j].pos[k],4)<=0.0) continue;
				for (m=0;m<4;m++) nav->peph[i].pos[k][m]=nav->peph[j].pos[k][m];
				for (m=0;m<4;m++) nav->peph[i].std[k][m]=nav->peph[j].std[k][m];
				for (m=0;m<4;m++) nav->peph[i].vel[k][m]=nav->peph[j].vel[k][m];
				for (m=0;m<4;m++) nav->peph[i].vst[k][m]=nav->peph[j].vst[k][m];
			}
		}
		else if (++i<j) nav->peph[i]=nav->peph[j];
	}
	nav->ne=i+1;

	trace(4,"combpeph: ne=%d\n",nav->ne);
}
/* read sp3 precise ephemeris file ---------------------------------------------
* read sp3 precise ephemeris/clock files and set them to navigation data
* args   : char   *file       I   sp3-c precise ephemeris file
*                                 (wind-card * is expanded)
*          nav_t  *nav        IO  navigation data
*          int    opt         I   options (1: only observed + 2: only predicted +
*                                 4: not combined)
* return : none
* notes  : see ref [1]
*          precise ephemeris is appended and combined
*          nav->peph and nav->ne must by properly initialized before calling the
*          function
*          only files with extensions of .sp3, .SP3, .eph* and .EPH* are read
*-----------------------------------------------------------------------------*/
extern void readsp3(const char *file, nav_t *nav, int opt)
{
	FILE *fp;
	gtime_t time={0};
	double bfact[2]={0};
	int i,j,n,ns,sats[MAXSAT]={0};
	char *efiles[MAXEXFILE],*ext,type=' ',tsys[4]="";

	trace(3,"readpephs: file=%s\n",file);

	for (i=0;i<MAXEXFILE;i++) {
		if (!(efiles[i]=(char *)malloc(1024))) {
			for (i--;i>=0;i--) free(efiles[i]);
			return;
		}
	}
	/* expand wild card in file path */
	n=expath(file,efiles,MAXEXFILE);

	for (i=j=0;i<n;i++) {
		if (!(ext=strrchr(efiles[i],'.'))) continue;

		if (!strstr(ext+1,"sp3")&&!strstr(ext+1,".SP3")&&
			!strstr(ext+1,"eph")&&!strstr(ext+1,".EPH")) continue;

		if (!(fp=fopen(efiles[i],"r"))) {
			trace(2,"sp3 file open error %s\n",efiles[i]);
			continue;
		}
		/* read sp3 header */
		ns=readsp3h(fp,&time,&type,sats,bfact,tsys);

		/* read sp3 body */
		readsp3b(fp,type,sats,ns,bfact,tsys,j++,opt,nav);

		fclose(fp);
	}
	for (i=0;i<MAXEXFILE;i++) free(efiles[i]);

	/* combine precise ephemeris */
	if (nav->ne>0) combpeph(nav,opt);
}
/* read satellite antenna parameters -------------------------------------------
* read satellite antenna parameters
* args   : char   *file       I   antenna parameter file
*          gtime_t time       I   time
*          nav_t  *nav        IO  navigation data
* return : status (1:ok,0:error)
* notes  : only support antex format for the antenna parameter file
*-----------------------------------------------------------------------------*/
extern int readsap(const char *file, gtime_t time, nav_t *nav)
{
	pcvs_t pcvs={0};
	pcv_t pcv0={0},*pcv;
	int i;

	trace(3,"readsap : file=%s time=%s\n",file,time_str(time,0));

	if (!readpcv(file,&pcvs)) return 0;

	for (i=0;i<MAXSAT;i++) {
		pcv=searchpcv(i+1,"",time,&pcvs);
		nav->pcvs[i]=pcv?*pcv:pcv0;
	}
	free(pcvs.pcv);
	return 1;
}
/* read dcb parameters file --------------------------------------------------*/
static int readdcbf(const char *file, nav_t *nav, const sta_t *sta)
{
	FILE *fp;
	double cbias;
	char buff[256],str1[32],str2[32]="";
	int i,j,sat,type=0;

	trace(3,"readdcbf: file=%s\n",file);

	if (!(fp=fopen(file,"r"))) {
		trace(2,"dcb parameters file open error: %s\n",file);
		return 0;
	}
	while (fgets(buff,sizeof(buff),fp)) {

		if      (strstr(buff,"DIFFERENTIAL (P1-P2) CODE BIASES")) type=1;
		else if (strstr(buff,"DIFFERENTIAL (P1-C1) CODE BIASES")) type=2;
		else if (strstr(buff,"DIFFERENTIAL (P2-C2) CODE BIASES")) type=3;

		if (!type||sscanf(buff,"%s %s",str1,str2)<1) continue;

		if ((cbias=str2num(buff,26,9))==0.0) continue;

		if (sta&&(!strcmp(str1,"G")||!strcmp(str1,"R"))) { /* receiver dcb */
			for (i=0;i<MAXRCV;i++) {
				if (!strcmp(sta[i].name,str2)) break;
			}
			if (i<MAXRCV) {
				j=!strcmp(str1,"G")?0:1;
				nav->rbias[i][j][type-1]=cbias*1E-9*CLIGHT; /* ns -> m */
			}
		}
		else if ((sat=satid2no(str1))) { /* satellite dcb */
			nav->cbias[sat-1][type-1]=cbias*1E-9*CLIGHT; /* ns -> m */
		}
	}
	fclose(fp);

	return 1;
}
/* read dcb parameters ---------------------------------------------------------
* read differential code bias (dcb) parameters
* args   : char   *file       I   dcb parameters file (wild-card * expanded)
*          nav_t  *nav        IO  navigation data
*          sta_t  *sta        I   station info data to inport receiver dcb
*                                 (NULL: no use)
* return : status (1:ok,0:error)
* notes  : currently only p1-c1 bias of code *.dcb file
*-----------------------------------------------------------------------------*/
extern int readdcb(const char *file, nav_t *nav, const sta_t *sta)
{
	int i,j,n;
	char *efiles[MAXEXFILE]={0};

	trace(3,"readdcb : file=%s\n",file);

	for (i=0;i<MAXSAT;i++) for (j=0;j<3;j++) {
		nav->cbias[i][j]=0.0;
	}
	for (i=0;i<MAXEXFILE;i++) {
		if (!(efiles[i]=(char *)malloc(1024))) {
			for (i--;i>=0;i--) free(efiles[i]);
			return 0;
		}
	}
	n=expath(file,efiles,MAXEXFILE);

	for (i=0;i<n;i++) {
		readdcbf(efiles[i],nav,sta);
	}
	for (i=0;i<MAXEXFILE;i++) free(efiles[i]);

	return 1;
}
/* add satellite fcb ---------------------------------------------------------*/
static int addfcb(nav_t *nav, gtime_t ts, gtime_t te, int sat,
				  const double *bias, const double *std)
{
	fcbd_t *nav_fcb;
	int i,j;

	if (nav->nf>0&&fabs(timediff(ts,nav->fcb[nav->nf-1].ts))<=1e-3) {
		for (i=0;i<3;i++) {
			nav->fcb[nav->nf-1].bias[sat-1][i]=bias[i];
			nav->fcb[nav->nf-1].std [sat-1][i]=std [i];
		}
		return 1;
	}
	if (nav->nf>=nav->nfmax) {
		nav->nfmax=nav->nfmax<=0?2048:nav->nfmax*2;
		if (!(nav_fcb=(fcbd_t *)realloc(nav->fcb,sizeof(fcbd_t)*nav->nfmax))) {
			free(nav->fcb); nav->nf=nav->nfmax=0;
			return 0;
		}
		nav->fcb=nav_fcb;
	}
	for (i=0;i<MAXSAT;i++) for (j=0;j<3;j++) {
		nav->fcb[nav->nf].bias[i][j]=nav->fcb[nav->nf].std[i][j]=0.0;
	}
	for (i=0;i<3;i++) {
		nav->fcb[nav->nf].bias[sat-1][i]=bias[i];
		nav->fcb[nav->nf].std [sat-1][i]=std [i];
	}
	nav->fcb[nav->nf  ].ts=ts;
	nav->fcb[nav->nf++].te=te;
	return 1;
}
/* read satellite fcb file ---------------------------------------------------*/
static int readfcbf(const char *file, nav_t *nav)
{
	FILE *fp;
	gtime_t ts,te;
	double ep1[6],ep2[6],bias[3]={0},std[3]={0};
	char buff[1024],str[32],*p;
	int sat;

	trace(3,"readfcbf: file=%s\n",file);

	if (!(fp=fopen(file,"r"))) {
		trace(2,"fcb parameters file open error: %s\n",file);
		return 0;
	}
	while (fgets(buff,sizeof(buff),fp)) {
		if ((p=strchr(buff,'#'))) *p='\0';
		if (sscanf(buff,"%lf/%lf/%lf %lf:%lf:%lf %lf/%lf/%lf %lf:%lf:%lf %s"
			"%lf %lf %lf %lf %lf %lf",ep1,ep1+1,ep1+2,ep1+3,ep1+4,ep1+5,
			ep2,ep2+1,ep2+2,ep2+3,ep2+4,ep2+5,str,bias,std,bias+1,std+1,
			bias+2,std+2)<17) continue;
		if (!(sat=satid2no(str))) continue;
		ts=epoch2time(ep1);
		te=epoch2time(ep2);
		if (!addfcb(nav,ts,te,sat,bias,std)) return 0;
	}
	fclose(fp);
	return 1;
}
/* compare satellite fcb -----------------------------------------------------*/
static int cmpfcb(const void *p1, const void *p2)
{
	fcbd_t *q1=(fcbd_t *)p1,*q2=(fcbd_t *)p2;
	double tt=timediff(q1->ts,q2->ts);
	return tt<-1E-3?-1:(tt>1E-3?1:0);
}
/* read satellite fcb data -----------------------------------------------------
* read satellite fractional cycle bias (dcb) parameters
* args   : char   *file       I   fcb parameters file (wild-card * expanded)
*          nav_t  *nav        IO  navigation data
* return : status (1:ok,0:error)
* notes  : fcb data appended to navigation data
*-----------------------------------------------------------------------------*/
extern int readfcb(const char *file, nav_t *nav)
{
	char *efiles[MAXEXFILE]={0};
	int i,n;

	trace(3,"readfcb : file=%s\n",file);

	for (i=0;i<MAXEXFILE;i++) {
		if (!(efiles[i]=(char *)malloc(1024))) {
			for (i--;i>=0;i--) free(efiles[i]);
			return 0;
		}
	}
	n=expath(file,efiles,MAXEXFILE);

	for (i=0;i<n;i++) {
		readfcbf(efiles[i],nav);
	}
	for (i=0;i<MAXEXFILE;i++) free(efiles[i]);

	if (nav->nf>1) {
		qsort(nav->fcb,nav->nf,sizeof(fcbd_t),cmpfcb);
	}
	return 1;
}
/* polynomial interpolation by Neville's algorithm ---------------------------*/
static double interppol(const double *x, double *y, int n)
{
	int i,j;

	for (j=1;j<n;j++) {
		for (i=0;i<n-j;i++) {
			y[i]=(x[i+j]*y[i]-x[i]*y[i+1])/(x[i+j]-x[i]);
		}
	}
	return y[0];
}
/* satellite position by precise ephemeris -----------------------------------*/
static int pephpos(gtime_t time, int sat, const nav_t *nav, double *rs,
				   double *dts, double *vare, double *varc)
{
	double t[NMAX+1],p[3][NMAX+1],c[2],*pos,std=0.0,s[3],sinl,cosl;
	int i,j,k,index;

	trace(4,"pephpos : time=%s sat=%2d\n",time_str(time,3),sat);

	rs[0]=rs[1]=rs[2]=dts[0]=0.0;

	if (nav->ne<NMAX+1||
		timediff(time,nav->peph[0].time)<-MAXDTE||
		timediff(time,nav->peph[nav->ne-1].time)>MAXDTE) {
			trace(3,"no prec ephem %s sat=%2d\n",time_str(time,0),sat);
			return 0;
	}
	/* binary search */
	for (i=0,j=nav->ne-1;i<j;) {
		k=(i+j)/2;
		if (timediff(nav->peph[k].time,time)<0.0) i=k+1; else j=k;
	}
	index=i<=0?0:i-1;

	/* polynomial interpolation for orbit */
	i=index-(NMAX+1)/2;
	if (i<0) i=0; else if (i+NMAX>=nav->ne) i=nav->ne-NMAX-1;

	for (j=0;j<=NMAX;j++) {
		t[j]=timediff(nav->peph[i+j].time,time);
		if (norm(nav->peph[i+j].pos[sat-1],3)<=0.0) {
			trace(3,"prec ephem outage %s sat=%2d\n",time_str(time,0),sat);
			return 0;
		}
	}
	for (j=0;j<=NMAX;j++) {
		pos=nav->peph[i+j].pos[sat-1];
#if 0
		p[0][j]=pos[0];
		p[1][j]=pos[1];
#else
		/* correciton for earh rotation ver.2.4.0 */
		sinl=sin(OMGE*t[j]);
		cosl=cos(OMGE*t[j]);
		p[0][j]=cosl*pos[0]-sinl*pos[1];
		p[1][j]=sinl*pos[0]+cosl*pos[1];
#endif
		p[2][j]=pos[2];
	}
	for (i=0;i<3;i++) {
		rs[i]=interppol(t,p[i],NMAX+1);
	}
	if (vare) {
		for (i=0;i<3;i++) s[i]=nav->peph[index].std[sat-1][i];
		std=norm(s,3);

		/* extrapolation error for orbit */
		if      (t[0   ]>0.0) std+=EXTERR_EPH*SQR(t[0   ])/2.0;
		else if (t[NMAX]<0.0) std+=EXTERR_EPH*SQR(t[NMAX])/2.0;
		*vare=SQR(std);
	}
	/* linear interpolation for clock */
	t[0]=timediff(time,nav->peph[index  ].time);
	t[1]=timediff(time,nav->peph[index+1].time);
	c[0]=nav->peph[index  ].pos[sat-1][3];
	c[1]=nav->peph[index+1].pos[sat-1][3];

	if (t[0]<=0.0) {
		if ((dts[0]=c[0])!=0.0) {
			std=nav->peph[index].std[sat-1][3]*CLIGHT-EXTERR_CLK*t[0];
		}
	}
	else if (t[1]>=0.0) {
		if ((dts[0]=c[1])!=0.0) {
			std=nav->peph[index+1].std[sat-1][3]*CLIGHT+EXTERR_CLK*t[1];
		}
	}
	else if (c[0]!=0.0&&c[1]!=0.0) {
		dts[0]=(c[1]*t[0]-c[0]*t[1])/(t[0]-t[1]);
		i=t[0]<-t[1]?0:1;
		std=nav->peph[index+i].std[sat-1][3]+EXTERR_CLK*fabs(t[i]);
	}
	else {
		dts[0]=0.0;
	}
	if (varc) *varc=SQR(std);
	return 1;
}
/* satellite clock by precise clock ------------------------------------------*/
static int pephclk(gtime_t time, int sat, const nav_t *nav, double *dts,
				   double *varc)
{
	double t[2],c[2],std;
	int i,j,k,index;

	trace(4,"pephclk : time=%s sat=%2d\n",time_str(time,3),sat);

	if (nav->nc<2||
		timediff(time,nav->pclk[0].time)<-MAXDTE||
		timediff(time,nav->pclk[nav->nc-1].time)>MAXDTE) {
			trace(3,"no prec clock %s sat=%2d\n",time_str(time,0),sat);
			return 1;
	}
	/* binary search */
	for (i=0,j=nav->nc-1;i<j;) {
		k=(i+j)/2;
		if (timediff(nav->pclk[k].time,time)<0.0) i=k+1; else j=k;
	}
	index=i<=0?0:i-1;

	/* linear interpolation for clock */
	t[0]=timediff(time,nav->pclk[index  ].time);
	t[1]=timediff(time,nav->pclk[index+1].time);
	c[0]=nav->pclk[index  ].clk[sat-1][0];
	c[1]=nav->pclk[index+1].clk[sat-1][0];

	if (t[0]<=0.0) {
		if ((dts[0]=c[0])==0.0) return 0;
		std=nav->pclk[index].std[sat-1][0]*CLIGHT-EXTERR_CLK*t[0];
	}
	else if (t[1]>=0.0) {
		if ((dts[0]=c[1])==0.0) return 0;
		std=nav->pclk[index+1].std[sat-1][0]*CLIGHT+EXTERR_CLK*t[1];
	}
	else if (c[0]!=0.0&&c[1]!=0.0) {
		dts[0]=(c[1]*t[0]-c[0]*t[1])/(t[0]-t[1]);
		i=t[0]<-t[1]?0:1;
		std=nav->pclk[index+i].std[sat-1][0]*CLIGHT+EXTERR_CLK*fabs(t[i]);
	}
	else {
		trace(3,"prec clock outage %s sat=%2d\n",time_str(time,0),sat);
		return 0;
	}
	if (varc) *varc=SQR(std);
	return 1;
}
/* satellite antenna phase center offset ---------------------------------------
* compute satellite antenna phase center offset in ecef
* args   : gtime_t time       I   time (gpst)
*          double *rs         I   satellite position and velocity (ecef)
*                                 {x,y,z,vx,vy,vz} (m|m/s)
*          int    sat         I   satellite number
*          nav_t  *nav        I   navigation data
*          double *dant       I   satellite antenna phase center offset (ecef)
*                                 {dx,dy,dz} (m) (iono-free LC value)
* return : none
*-----------------------------------------------------------------------------*/
extern void satantoff(gtime_t time, const double *rs, int sat, const nav_t *nav,
					  double *dant)
{
	const double *lam=nav->lam[sat-1];
	const pcv_t *pcv=nav->pcvs+sat-1;
	double ex[3],ey[3],ez[3],es[3],r[3],rsun[3],gmst,erpv[5]={0};
	double gamma,C1,C2,dant1,dant2;
	int i,j=0,k=1;

	trace(4,"satantoff: time=%s sat=%2d\n",time_str(time,3),sat);

	/* sun position in ecef */
	sunmoonpos(gpst2utc(time),erpv,rsun,NULL,&gmst);

	/* unit vectors of satellite fixed coordinates */
	for (i=0;i<3;i++) r[i]=-rs[i];
	if (!normv3(r,ez)) return;
	for (i=0;i<3;i++) r[i]=rsun[i]-rs[i];
	if (!normv3(r,es)) return;
	cross3(ez,es,r);
	if (!normv3(r,ey)) return;
	cross3(ey,ez,ex);

	if (NFREQ>=3&&(satsys(sat,NULL)&(SYS_GAL|SYS_SBS))) k=2;

	if (NFREQ<2||lam[j]==0.0||lam[k]==0.0) return;

	gamma=SQR(lam[k])/SQR(lam[j]);
	C1=gamma/(gamma-1.0);
	C2=-1.0 /(gamma-1.0);

	/* iono-free LC */
	for (i=0;i<3;i++) {
		dant1=pcv->off[j][0]*ex[i]+pcv->off[j][1]*ey[i]+pcv->off[j][2]*ez[i];
		dant2=pcv->off[k][0]*ex[i]+pcv->off[k][1]*ey[i]+pcv->off[k][2]*ez[i];
		dant[i]=C1*dant1+C2*dant2;
	}
}
/* satellite position/clock by precise ephemeris/clock -------------------------
* compute satellite position/clock with precise ephemeris/clock
* args   : gtime_t time       I   time (gpst)
*          int    sat         I   satellite number
*          nav_t  *nav        I   navigation data
*          int    opt         I   sat postion option
*                                 (0: center of mass, 1: antenna phase center)
*          double *rs         O   sat position and velocity (ecef)
*                                 {x,y,z,vx,vy,vz} (m|m/s)
*          double *dts        O   sat clock {bias,drift} (s|s/s)
*          double *var        IO  sat position and clock error variance (m)
*                                 (NULL: no output)
* return : status (1:ok,0:error or data outage)
* notes  : clock includes relativistic correction but does not contain code bias
*          before calling the function, nav->peph, nav->ne, nav->pclk and
*          nav->nc must be set by calling readsp3(), readrnx() or readrnxt()
*          if precise clocks are not set, clocks in sp3 are used instead
*-----------------------------------------------------------------------------*/
extern int peph2pos(gtime_t time, int sat, const nav_t *nav, int opt,
					double *rs, double *dts, double *var)
{
	double rss[3],rst[3],dtss[1],dtst[1],dant[3]={0},vare=0.0,varc=0.0,tt=1E-3;
	int i;

	trace(4,"peph2pos: time=%s sat=%2d opt=%d\n",time_str(time,3),sat,opt);

	if (sat<=0||MAXSAT<sat) return 0;

	/* satellite position and clock bias */
	if (!pephpos(time,sat,nav,rss,dtss,&vare,&varc)||
		!pephclk(time,sat,nav,dtss,&varc)) return 0;

	time=timeadd(time,tt);
	if (!pephpos(time,sat,nav,rst,dtst,NULL,NULL)||
		!pephclk(time,sat,nav,dtst,NULL)) return 0;

	/* satellite antenna offset correction */
	if (opt) {
		satantoff(time,rss,sat,nav,dant);
	}
	for (i=0;i<3;i++) {
		rs[i  ]=rss[i]+dant[i];
		rs[i+3]=(rst[i]-rss[i])/tt;
	}
	/* relativistic effect correction */
	if (dtss[0]!=0.0) {
		dts[0]=dtss[0]-2.0*dot(rs,rs+3,3)/CLIGHT/CLIGHT;
		dts[1]=(dtst[0]-dtss[0])/tt;
	}
	else { /* no precise clock */
		dts[0]=dts[1]=0.0;
	}
	if (var) *var=vare+varc;

	return 1;
}
