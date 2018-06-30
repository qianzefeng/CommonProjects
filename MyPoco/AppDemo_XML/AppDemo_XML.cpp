#include "DOMParser.hpp"
#include "DOMWriter.hpp"
#include "PrettyPrint.hpp"
#include "SAXParser.hpp"

int main(int argc, char** argv)
{
	DOMParser_main(argc, argv);
	DOMWriter_main(argc, argv);
	PrettyPrint_main(argc, argv);
	SAXParser_main(argc, argv);

	system("pause");
	return 0;
}