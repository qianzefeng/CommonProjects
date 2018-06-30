//
// Handler.cpp
//
// Library: JSON
// Package: JSON
// Module:  Handler
//
// Copyright (c) 2012, Applied Informatics Software Engineering GmbH.
// and Contributors.
//
// SPDX-License-Identifier:	BSL-1.0
//


#include "JSON/Handler.h"
#include "JSON/Object.h"


namespace Poco {
namespace JSON {


Handler::Handler()
{
}


Handler::~Handler()
{
}


Dynamic::Var Handler::asVar() const
{
	return Dynamic::Var();
}


Poco::DynamicStruct Handler::asStruct() const
{
	return Poco::DynamicStruct();
}


} } // namespace Poco::JSON
