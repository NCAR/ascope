/*
 * ScopeTSSource.cpp
 *
 *  Created on: Jan 28, 2009
 *      Author: martinc
 */

#include <iostream>

#include "ScopeTSReader.h"

//////////////////////////////////////////////////////////////////////////////
ScopeTSReader::ScopeTSReader(DDSSubscriber& subscriber, std::string topicName, double outputRate):
	QtTSReader(subscriber, topicName, outputRate)
{
	connect(this, SIGNAL(newItem(ProfilerDDS::TimeSeries*)),
			this, SLOT(newItemSlot(ProfilerDDS::TimeSeries*)));
}

//////////////////////////////////////////////////////////////////////////////
ScopeTSReader::~ScopeTSReader() {
}

//////////////////////////////////////////////////////////////////////////////
void
ScopeTSReader::newItemSlot(ProfilerDDS::TimeSeries* pItem) {
	std::cout << "returning a time series" << std::endl;
	returnItem(pItem);
}
