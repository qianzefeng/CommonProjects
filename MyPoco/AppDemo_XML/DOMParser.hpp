//
// DOMParser.cpp
//
// This sample demonstrates the DOMParser, AutoPtr and
// NodeIterator classes.
//
// Copyright (c) 2004-2006, Applied Informatics Software Engineering GmbH.
// and Contributors.
//
// SPDX-License-Identifier:	BSL-1.0
//


#include "XML/DOM/DOMParser.h"
#include "XML/DOM/Document.h"
#include "XML/DOM/NodeIterator.h"
#include "XML/DOM/NodeFilter.h"
#include "XML/DOM/AutoPtr.h"
#include "XML/SAX/InputSource.h"
#include "Poco/Exception.h"
#include <iostream>


using Poco::XML::DOMParser;
using Poco::XML::InputSource;
using Poco::XML::Document;
using Poco::XML::NodeIterator;
using Poco::XML::NodeFilter;
using Poco::XML::Node;
using Poco::XML::AutoPtr;
using Poco::Exception;


int DOMParser_main(int argc, char** argv)
{
	// Parse an XML document from standard input
	// and use a NodeIterator to print out all nodes.
	
	InputSource src(std::cin);
	try
	{
		DOMParser parser;
		AutoPtr<Document> pDoc = parser.parse(&src);
		
		NodeIterator it(pDoc, NodeFilter::SHOW_ALL);
		Node* pNode = it.nextNode();
		while (pNode)
		{
			std::cout << pNode->nodeName() << ":" << pNode->nodeValue() << std::endl;
			pNode = it.nextNode();
		}
	}
	catch (Exception& exc)
	{
		std::cerr << exc.displayText() << std::endl;
	}
	
	return 0;
}
