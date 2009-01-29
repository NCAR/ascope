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
	QWidget(parent)
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

	std::cout << "Returning a time series."
	    << " tsdata size:" << size
	    << " Current dropped samples:" << _tsReader->droppedSamples()
		<< " Current number of samples:" << _tsReader->numSamples()
		<< std::endl;

	_tsReader->returnItem(pItem);
}
