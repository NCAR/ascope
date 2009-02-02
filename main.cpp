#include <QApplication>
#include "ProfilerScope.h"


int
main (int argc, char** argv) {
	QApplication app(argc, argv);

    // create a dialog to serve as parent
    QDialog* dialog = new QDialog;
	dialog->show();

	// create the scope
	ProfilerScope scope(dialog);

	return app.exec();
}
