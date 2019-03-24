#ifndef  _RTKLIB_APP_H_
#define  _RTKLIB_APP_H_

#ifdef __cplusplus
extern "C" {
#endif


/*
atest1:
./rnx2rtcm $(OPT1) -typ 1002,1019,1033 -sta 111 -out out/rtcm_1002.rtcm3
test2:
./rnx2rtcm $(OPT1) -typ 1004,1019,1033 -sta 222 -out out/rtcm_1004.rtcm3
test3:
./rnx2rtcm $(OPT1) -typ 1010,1020,1033 -sta 333 -out out/rtcm_1010.rtcm3
test4:
./rnx2rtcm $(OPT1) -typ 1012,1020,1033 -sta 444 -out out/rtcm_1012.rtcm3
test5:
./rnx2rtcm $(OPT1) -typ 1044 -sta 555 -out out/rtcm_1044.rtcm3
test6:
./rnx2rtcm $(OPT1) -typ 1074 -sta 666 -out out/rtcm_1074.rtcm3
test7:
./rnx2rtcm $(OPT2) -typ 1074,1084,1094,1104,1114,1019,1020,1044,1045,1046 -sta 444 -out out/rtcm_1074_jav.rtcm3
test8:
./rnx2rtcm $(OPT2) -typ 1075,1085,1095,1105,1115,1019,1020,1044,1045,1046 -sta 555 -out out/rtcm_1075_jav.rtcm3
test9:
./rnx2rtcm $(OPT2) -typ 1076,1086,1096,1106,1116,1019,1020,1044,1045,1046 -sta 666 -out out/rtcm_1076_jav.rtcm3
test10:
./rnx2rtcm $(OPT2) -typ 1077,1087,1097,1107,1117,1019,1020,1044,1045,1046 -sta 777 -out out/rtcm_1077_jav.rtcm3
test11:
./rnx2rtcm $(OPT3) -typ 1077,1087,1019,1020 -sta 888 -out out/rtcm_1077_gras.rtcm3
*/
extern int rnx2rtcm_main(int argc, char **argv);

extern int t_time_main(void);

extern int t_rinex_main(int argc, char **argv);

extern int t_stec_main(int argc, char **argv);

extern int t_preceph_main(int argc, char **argv);

extern int t_tle_main(int argc, char **argv);

extern int t_atmos_main(void);

extern int t_coord_main(void);

extern int t_geoid_main(void);

extern int t_ppp_main(void);

extern int t_misc_main(void);

extern int t_gloeph_main(int argc, char **argv);

extern int t_ionex_main(int argc, char **argv);

extern int t_lambda_main(void);

extern int t_matrix_main(void);


#ifdef __cplusplus
}
#endif
#endif /* RTKLIB_H */