/*
 * ProfilerScope.cpp
 *
 *  Created on: Jan 29, 2009
 *      Author: martinc
 */


#include "ProfilerScope.h"

/////////////////////////////////////////////////////////////////////////
ProfilerScope::ProfilerScope(TSReader* reader, QWidget* parent):
	_tsReader(reader),
	ScopePlot(parent)
{

}

/////////////////////////////////////////////////////////////////////////
ProfilerScope::~ProfilerScope() {

}

/////////////////////////////////////////////////////////////////////////
void ProfilerScope::newTSItemSlot(ProfilerDDS::TimeSeries* pItem) {
	int size = pItem->tsdata.length();
	int gates = pItem->hskp.gates;
	int channels = pItem->hskp.numChannels;
	int tsLength = pItem->hskp.tsLength;
	
	std::vector<double> I, Q;
	I.resize(gates*tsLength);
	Q.resize(gates*tsLength);
	
	int index = 0;
	for (unsigned int g = 0; g < gates; g++) {
		for (unsigned int t = 0; t < tsLength; t++) {
			I[g*tsLength + t] = pItem->tsdata[index++];
			Q[g*tsLength + t] = pItem->tsdata[index++];
		}
	}

	this->TimeSeries(I, Q, -100.0, 100.0, 1.0, "i", "I&Q");
	std::cout << "Returning a time series."
	    << " tsdata size:" << size
	    << " Current dropped samples:" << _tsReader->droppedSamples()
		<< " Current number of samples:" << _tsReader->numSamples()
		<< std::endl;

	_tsReader->returnItem(pItem);
}
