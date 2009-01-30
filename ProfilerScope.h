/*
 * ProfilerScope.h
 *
 *  Created on: Jan 29, 2009
 *      Author: martinc
 */

#ifndef PROFILERSCOPE_H_
#define PROFILERSCOPE_H_

#include <QWidget>
#include "ScopePlot.h"
#include "TSReader.h"

class ProfilerScope: public ScopePlot {
	Q_OBJECT
public:
	ProfilerScope(TSReader* reader, QWidget* parent=0);
	virtual ~ProfilerScope();

public slots:
	void newTSItemSlot(ProfilerDDS::TimeSeries*);

protected:
	TSReader* _tsReader;

};

#endif /* PROFILERSCOPE_H_ */
