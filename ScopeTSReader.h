/*
 * ScopeTSSource.h
 *
 *  Created on: Jan 28, 2009
 *      Author: martinc
 */

#ifndef SCOPETSREADER_H_
#define SCOPETSREADER_H_

#include "DDSSubscriber.h"
#include "QtTSReader.h"

class ScopeTSReader: public QtTSReader {
	Q_OBJECT
public:
	/// Constructor
	/// @param subscriber The DDS subscriber
	/// @param topicName The DDS topic name
	/// @param outputRate The desired maximum sample output rate, in Hz
	ScopeTSReader(DDSSubscriber& subscriber, std::string topicName, double outputRate);
	/// Destructor
	virtual ~ScopeTSReader();

public slots:
    /// This slot will be connected to the newItem()
    /// signal from QtTSReader().
    /// @todo It will do what with each sample? However,
    /// it will return the item.
    void newItemSlot(ProfilerDDS::TimeSeries* pItem);


};
#endif /* SCOPETSREADER_H_ */
