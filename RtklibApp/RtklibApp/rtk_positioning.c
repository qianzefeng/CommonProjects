#include "rtk_positioning.h"
#include "rtk_mat.h"
#include "rtk_coord.h"
#include "rtk_sat.h"
#include "rtk_trace.h"




/* satellite carrier wave length -----------------------------------------------
* get satellite carrier wave lengths
* args   : int    sat       I   satellite number
*          int    frq       I   frequency index (0:L1,1:L2,2:L5/3,...)
*          nav_t  *nav      I   navigation messages
* return : carrier wave length (m) (0.0: error)
*-----------------------------------------------------------------------------*/
extern double satwavelen(int sat, int frq, const nav_t *nav)
{
	const double freq_glo[]={FREQ1_GLO,FREQ2_GLO,FREQ3_GLO};
	const double dfrq_glo[]={DFRQ1_GLO,DFRQ2_GLO,0.0};
	int i,sys=satsys(sat, NULL);

	if (sys==SYS_GLO) 
	{
		if (0<=frq&&frq<=2) 
		{
			for (i=0;i<nav->ng;i++) 
			{
				if (nav->geph[i].sat!=sat) continue;
				return CLIGHT/(freq_glo[frq]+dfrq_glo[frq]*nav->geph[i].frq);
			}
		}
	}
	else if (sys==SYS_CMP) 
	{
		if      (frq==0) return CLIGHT/FREQ1_CMP; /* B1 */
		else if (frq==1) return CLIGHT/FREQ2_CMP; /* B3 */
		else if (frq==2) return CLIGHT/FREQ3_CMP; /* B2 */
	}
	else 
	{
		if      (frq==0) return CLIGHT/FREQ1; /* L1/E1 */
		else if (frq==1) return CLIGHT/FREQ2; /* L2 */
		else if (frq==2) return CLIGHT/FREQ5; /* L5/E5a */
		else if (frq==3) return CLIGHT/FREQ6; /* L6/LEX */
		else if (frq==4) return CLIGHT/FREQ7; /* E5b */
		else if (frq==5) return CLIGHT/FREQ8; /* E5a+b */
	}
	return 0.0;
}



/* satellite azimuth/elevation angle -------------------------------------------
* compute satellite azimuth/elevation angle
* args   : double *pos      I   geodetic position {lat,lon,h} (rad,m)
*          double *e        I   receiver-to-satellilte unit vevtor (ecef)
*          double *azel     IO  azimuth/elevation {az,el} (rad) (NULL: no output)
*                               (0.0<=azel[0]<2*pi,-pi/2<=azel[1]<=pi/2)
* return : elevation angle (rad)
*-----------------------------------------------------------------------------*/
extern double satazel(const double *pos, const double *e, double *azel)
{
	double az=0.0,el=PI/2.0,enu[3];

	if (pos[2]>-RE_WGS84) 
	{
		ecef2enu(pos,e,enu);
		az=dot(enu,enu,2)<1E-12?0.0:atan2(enu[0],enu[1]);
		if (az<0.0) az+=2*PI;
		el=asin(enu[2]);
	}
	if (azel) {azel[0]=az; azel[1]=el;}
	return el;
}

/* geometric distance ----------------------------------------------------------
* compute geometric distance and receiver-to-satellite unit vector
* args   : double *rs       I   satellilte position (ecef at transmission) (m)
*          double *rr       I   receiver position (ecef at reception) (m)
*          double *e        O   line-of-sight vector (ecef)
* return : geometric distance (m) (0>:error/no satellite position)
* notes  : distance includes sagnac effect correction
*-----------------------------------------------------------------------------*/
extern double geodist(const double *rs, const double *rr, double *e)
{
	double r;
	int i;

	if (norm(rs,3)<RE_WGS84) return -1.0;
	for (i=0;i<3;i++) e[i]=rs[i]-rr[i];
	r=norm(e,3);
	for (i=0;i<3;i++) e[i]/=r;
	return r+OMGE*(rs[0]*rr[1]-rs[1]*rr[0])/CLIGHT;
}

/* compute dops ----------------------------------------------------------------
* compute DOP (dilution of precision)
* args   : int    ns        I   number of satellites
*          double *azel     I   satellite azimuth/elevation angle (rad)
*          double elmin     I   elevation cutoff angle (rad)
*          double *dop      O   DOPs {GDOP,PDOP,HDOP,VDOP}
* return : none
* notes  : dop[0]-[3] return 0 in case of dop computation error
*-----------------------------------------------------------------------------*/
#define SQRT(x)     ((x)<0.0?0.0:sqrt(x))

extern void dops(int ns, const double *azel, double elmin, double *dop)
{
	double H[4*MAXSAT],Q[16],cosel,sinel;
	int i,n;

	for (i=0;i<4;i++) dop[i]=0.0;
	for (i=n=0;i<ns&&i<MAXSAT;i++) 
	{
		if (azel[1+i*2]<elmin||azel[1+i*2]<=0.0) continue;
		cosel=cos(azel[1+i*2]);
		sinel=sin(azel[1+i*2]);
		H[  4*n]=cosel*sin(azel[i*2]);
		H[1+4*n]=cosel*cos(azel[i*2]);
		H[2+4*n]=sinel;
		H[3+4*n++]=1.0;
	}
	if (n<4) return;

	matmul("NT",4,4,n,1.0,H,H,0.0,Q);
	if (!matinv(Q,4)) 
	{
		dop[0]=SQRT(Q[0]+Q[5]+Q[10]+Q[15]); /* GDOP */
		dop[1]=SQRT(Q[0]+Q[5]+Q[10]);       /* PDOP */
		dop[2]=SQRT(Q[0]+Q[5]);             /* HDOP */
		dop[3]=SQRT(Q[10]);                 /* VDOP */
	}
}

/* carrier smoothing -----------------------------------------------------------
* carrier smoothing by Hatch filter
* args   : obs_t  *obs      IO  raw observation data/smoothed observation data
*          int    ns        I   smoothing window size (epochs)
* return : none
*-----------------------------------------------------------------------------*/
extern void csmooth(obs_t *obs, int ns)
{
	double Ps[2][MAXSAT][NFREQ]={{{0}}},Lp[2][MAXSAT][NFREQ]={{{0}}},dcp;
	int i,j,s,r,n[2][MAXSAT][NFREQ]={{{0}}};
	obsd_t *p;

	trace(3,"csmooth: nobs=%d,ns=%d\n",obs->n,ns);

	for (i=0;i<obs->n;i++) 
	{
		p=&obs->data[i]; s=p->sat; r=p->rcv;
		for (j=0;j<NFREQ;j++) 
		{
			if (s<=0||MAXSAT<s||r<=0||2<r) continue;
			if (p->P[j]==0.0||p->L[j]==0.0) continue;
			if (p->LLI[j]) n[r-1][s-1][j]=0;
			if (n[r-1][s-1][j]==0) Ps[r-1][s-1][j]=p->P[j];
			else 
			{
				dcp=lam_carr[j]*(p->L[j]-Lp[r-1][s-1][j]);
				Ps[r-1][s-1][j]=p->P[j]/ns+(Ps[r-1][s-1][j]+dcp)*(ns-1)/ns;
			}
			if (++n[r-1][s-1][j]<ns) p->P[j]=0.0; else p->P[j]=Ps[r-1][s-1][j];
			Lp[r-1][s-1][j]=p->L[j];
		}
	}
}
