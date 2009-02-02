#include <QApplication>
#include "ProfilerScope.h"


int
main (int argc, char** argv) {
	QApplication app(argc, argv);

 	// create the scope
	ProfilerScope scope;
	scope.show();

	return app.exec();
}
