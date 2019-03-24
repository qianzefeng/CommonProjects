#include "RTKLIB_App.h"


int main(int argc, char **argv)
{
	t_time_main();
	t_rinex_main(argc, argv);
	rnx2rtcm_main(argc, argv);
	
	system("pause");
	return 0;
}