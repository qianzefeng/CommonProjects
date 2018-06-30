#include "ActiveMethod.hpp"
#include "Activity.hpp"
#include "base64decode.hpp"
#include "base64encode.hpp"
#include "BinaryReaderWriter.hpp"
#include "DateTime.hpp"
#include "deflate.hpp"
#include "grep.hpp"
#include "hmacmd5.hpp"
#include "inflate.hpp"
#include "LineEndingConverter.hpp"
#include "Logger.hpp"
#include "LogRotation.hpp"
#include "md5.hpp"
#include "NotificationQueue.hpp"
#include "StringTokenizer.hpp"
#include "URI.hpp"
#include "uuidgen.hpp"
#include "Timer.hpp"



int main(int argc, char** argv)
{
	ActiveMethod_Test(argc, argv);
	Activity_Test(argc, argv);
	base64decode_Test(argc, argv);
	base64encode_Test(argc, argv);
	BinaryReaderWriter_Test(argc, argv);
	DateTime_Test(argc, argv);
	deflate_Test(argc, argv);
	grep_Test(argc, argv);
	hmacmd5_Test(argc, argv);
	inflate_Test(argc, argv);
	LineEndingConverter_Test(argc, argv);
	Logger_Test(argc, argv);
	LogRotation_Test(argc, argv);
	md5_Test(argc, argv);
	NotificationQueue_Test(argc, argv);
	StringTokenizer_Test(argc, argv);
	URI_Test(argc, argv);
	uuidgen_Test(argc, argv);
	Timer_Test(argc, argv);

	system("pause");
	return 0;
}