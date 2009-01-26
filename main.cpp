/*
 * main.cpp
 *
 *  Created on: Jan 21, 2009
 *      Author: martinc
 */

//#include "DDSQtReader.h"

#include "DDSSubscriber.h"
#include "QtTSReader.h"

int
main (int argc, char** argv) {

	DDSSubscriber subscriber(argc, argv);

	QtTSReader test(subscriber, std::string("ProfilerTS"), 100.0);

	return 0;
}
