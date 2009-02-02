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

/// ScopeTSReader subclasses QtTSReader and monitors the
/// incoming time series sample stream, emitting diagnostic
/// signals.
///
/// The real end user of the time series samples must connect
/// to the QtTSReader::newItem(ProfilerDDS::TimeSeries*) signal
/// and harvest the DDS samples, and then return them via returnItem()
/// when finished.
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

signals:

private slots:
    /// This slot will be connected to the newItem()
    /// signal from QtTSReader(). It provides a place
    /// place for ScopeTSReader to monitor or react to
    /// the time series sample stream. However, it is
    /// not repsonsible for returning the pItem.
    void newItemSlot(ProfilerDDS::TimeSeries* pItem);
    /// Use this slot to return an item
    /// @param pItem the item to be returned.
    void returnItemSlot(ProfilerDDS::TimeSeries* pItem);
    

};
#endif /* SCOPETSREADER_H_ */
