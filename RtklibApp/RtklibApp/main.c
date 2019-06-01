#include <assert.h>

#include "rtk_time.h"
#include "rtk_rtcm.h"
#include "rtk_rinex.h"
#include "rtk_io.h"
#include "tst_rnx2rtcm.h"

void rtk_time_fun1();
void rtk_time_fun2();
void rtk_time_fun3();
void rtk_decode_rtcm3();
void rtk_encode_rtcm3();


void GetModulePath(char strDir[260])
{
	char strPath[260] = {'\0'};
	if (GetModuleFileNameA(NULL, strPath, 259) > 0)
	{
		char* strFile = NULL;
		if( NULL != strrchr(strPath, '\\') )
		{
			strFile = strrchr(strPath, '\\' ) + 1;
			memcpy(strDir, strPath, strlen(strPath) - strlen(strFile));
		}
	}
}

void main(int argc, char **argv)
{
	char strDir[260] = {'\0'};
	GetModulePath(strDir);

	rtk_time_fun1();
	rtk_time_fun2();
	rtk_time_fun3();

	//rtk_decode_rtcm3();
	rtk_encode_rtcm3();

	getchar();
}


void rtk_encode_rtcm3()
{
	gtime_t ts={0},te={0};
	obs_t  obs={0};
	nav_t  nav={0};
	sta_t  sta={{0}};
	double es[6]={0}, ee[6]={0}, tint=0.0;
	char *infile[16]={0};
	char *outfile="rtklib_encode_rtcm.dat";
	int i, m = 0, type[16], staid=0, ret=0;

	type[0] = 1074;
	type[1] = 1084;
	type[2] = 1019;
	type[3] = 1020;
	type[4] = 1006;
	type[5] = 1033;
	m = 6;

	infile[0] = "522816178.15O";
	infile[1] = "522816178.15N";
	infile[2] = "522816178.15G";
	infile[3] = "522816178.15R";

	/* read rinex files */
	for (i = 0; i < m; i++) 
	{
		readrnxt(infile[i], 0, ts, te, tint, "", &obs, &nav, &sta);
	}

	sortobs(&obs);
	uniqnav(&nav);

	/* convert to rtcm messages */
	if (!conv_rtcm(type, m, outfile, &obs, &nav, &sta, staid)) ret=-1;

	free(obs.data);
	freenav(&nav, 0xFF);
}

void rtk_decode_rtcm3()
{
	rtcm_t rtcm = {0};
	FILE * pFile = NULL;
	char strPath[256] = {'\0'};
	printf("Please input rtcm file path:\n");
	scanf("%s", strPath);
	pFile = fopen(strPath, "rb");
	if (NULL != pFile)
	{
		init_rtcm(&rtcm);
		while (1)
		{
			if(input_rtcm3f(&rtcm, pFile) < -1)
			{
				break;
			}
		}
	}
	else
	{
		printf("Open file fail!\n");
	}
}

void rtk_time_fun3()
{
	gtime_t t;
	double  ep[6];
	int     w, week;
	double  time,tt;
	t = gpst2time(0, 0.0); 
	time2epoch(t, ep);
	assert(ep[0] == 1980 && ep[1] == 1 && ep[2] == 6&& ep[3] == 0 && ep[4] == 0 && ep[5] == 0.0);
	t = gpst2time(1400, 86400.0); 
	time2epoch(t, ep);
	assert(ep[0] == 2006 && ep[1] == 11 && ep[2] == 6 && ep[3] == 0 && ep[4] == 0 && ep[5] == 0.0);
	t = gpst2time(1400, 86400.0*7-1.0); 
	time2epoch(t,ep);
	assert(ep[0] == 2006 && ep[1] == 11 && ep[2] == 11 && ep[3] == 23 && ep[4] == 59 && ep[5] == 59.0);
	t = gpst2time(1400, 86400.0*7); 
	time2epoch(t,ep);
	assert(ep[0] == 2006 && ep[1] == 11&& ep[2] == 12 && ep[3] == 0 && ep[4] == 0 && ep[5] == 0.0);
	t = gpst2time(1401,0.0); 
	time2epoch(t,ep);
	assert(ep[0] == 2006 && ep[1] == 11 && ep[2] == 12 && ep[3] == 0 && ep[4] == 0 && ep[5] == 0.0);
	for (w = 1000; w < 1100; w++)
	{
		for(time = 0.0; time < 86400.0 * 7; time += 3600.0)
		{
			t  = gpst2time(w, time); 
			tt = time2gpst(t, &week);
			assert(tt == time && week == w);
		}
	}
	printf("%s rtk_time_fun3 : OK\n", __FILE__);
}

void rtk_time_fun2()
{
	double eps[]={1970, 1, 1, 0, 0, 0.000000};
	double ep0[]={1980, 1, 6, 0, 0, 0.000000};
	double ep1[]={2004, 2,28, 2, 0,59.999999};
	double ep2[]={2037,10, 1, 0, 0, 0.000000};
	int year, month, day, mday[]={31,28,31,30,31,30,31,31,30,31,30,31};
	gtime_t t;
	double ep[6];

	t = epoch2time(eps);
	time2epoch(t,ep);
	assert(ep[0] == 1970 && ep[1] == 1 && ep[2] == 1 && ep[3] == 0 && ep[4] == 0 && ep[5] == 0.0);
	t = epoch2time(ep0); 
	time2epoch(t,ep);
	assert(ep[0] == 1980 && ep[1] == 1 && ep[2] == 6 && ep[3] == 0 && ep[4] == 0 && ep[5] == 0.0);
	t = epoch2time(ep1);
	time2epoch(t,ep);
	assert(ep[0] == 2004 && ep[1] == 2 && ep[2] == 28 && ep[3] == 2 && ep[4] == 0 && fabs(ep[5]-59.999999) < 1E-14);
	t = epoch2time(ep2);
	time2epoch(t,ep);
	assert(ep[0] == 2037 && ep[1] == 10 && ep[2] == 1 && ep[3] == 0 && ep[4] == 0 && ep[5] == 0.0);

	for (year = 1970; year <= 2037; year++) 
	{
		mday[1] = year % 4 == 0 ? 29:28;
		for (month = 1; month <= 12; month++) 
		{
			for (day = 1; day <= mday[month-1]; day++) 
			{
				if (year == 1970 && month == 1 && day == 1)
				{
					continue;
				}
				ep0[0] = year; 
				ep0[1] = month;
				ep0[2] = day;
				t = epoch2time(ep0); 
				time2epoch(t, ep);
				assert(ep[0] == ep0[0] && ep[1] == ep0[1]&& ep[2] == ep0[2]);
				assert(ep[3] == 0.0    && ep[4] == 0.0   && ep[5] == 0.0);
			}
		}
	}
	printf("%s rtk_time_fun2 : OK\n",__FILE__);
}

void rtk_time_fun1()
{
	double  a;
	int     s;
	gtime_t t;
	double ep[6];
	char   *p = "$$3456.789................";
	char   str[80]= "....3D45......................";
	a = str2num(p, 2, 7);
	assert(fabs(a - 3456.78 ) < 1E-12);
	a = str2num(str, 4, 4);
	assert(fabs(a - 3E45    ) < 1E-30);


	strcpy(str, "....2004 1 1 0 1 2.345........");
	s = str2time(str, 0, 0, &t);//Ê§°Ü
	assert(s < 0);
	s = str2time(str, 4,17, &t);
	time2epoch(t, ep);
	assert(!s && ep[0] == 2004 && ep[1] == 1 && ep[2] == 1 && ep[3] == 0 && ep[4] == 1 && fabs(ep[5]-2.34) < 1E-15);

	strcpy(str, "....  37 12 31 1 2 3 .........");
	s = str2time(str, 4, 20, &t); 
	time2epoch(t, ep);
	assert(!s && ep[0] == 2037 && ep[1] == 12 && ep[2] == 31 && ep[3] == 1 && ep[4] == 2 && ep[5] == 3);

	printf("%s rtk_time_fun1 : OK\n",__FILE__);
}