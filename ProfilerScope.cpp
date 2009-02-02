#include "ProfilerScope.h"
#include "ScopePlot.h"
#include "Knob.h"
//#include "SvnVersion.h"

#include <QMessageBox>
#include <QButtonGroup>
#include <QLabel>
#include <QTimer>
#include <QSpinBox>
#include <QLCDNumber>
#include <QSlider>
#include <QLayout>
#include <QTabWidget>
#include <QWidget>
#include <QRadioButton>
#include <QButtonGroup>
#include <QFrame>
#include <QPushButton>
#include <QPalette>
#include <QDateTime>
#include <QFileDialog>
#include <string>
#include <algorithm>

#include <iostream>
#include <time.h>

#include <qwt_wheel.h>

//////////////////////////////////////////////////////////////////////
ProfilerScope::ProfilerScope(
        QWidget* parent) :
    QWidget(parent),
    _statsUpdateInterval(5),
            _timeSeriesPlot(TRUE), _config("NCAR", "ProfilerScope"),
            _paused(false), _zeroMoment(0.0){
    // Set up our form
    setupUi(this);

    return;

    // Initialize fft calculations
    initFFT();

    // get our title from the coniguration
    std::string title = _config.getString("title", "ProfilerScope");
    title += " ";
    //title += SvnVersion::revision();
    parent->setWindowTitle(title.c_str());

    // initialize running statistics
    for (int i = 0; i < 3; i++) {
        //		_pulseCount[i]	= 0;
        //		_prevPulseCount[i] = 0;
        _errorCount[i] = 0;
        _lastPulseNum[i] = 0;
    }

    // configure the channel selections
    _channel = 1;
    _chan1->setChecked(true);
    _chan2->setChecked(false);
    _chan3->setChecked(false);
    _chan4->setChecked(false);
    QButtonGroup* channelButtonGroup = new QButtonGroup();
    channelButtonGroup->addButton(_chan1, 1);
    channelButtonGroup->addButton(_chan2, 2);
    channelButtonGroup->addButton(_chan3, 3);
    channelButtonGroup->addButton(_chan4, 4);

    // connect the controls
    connect(_autoScale, SIGNAL(released()), this, SLOT(autoScaleSlot()));
    connect(_gainKnob, SIGNAL(valueChanged(double)), this, SLOT(gainChangeSlot(double)));
    connect(_up, SIGNAL(released()), this, SLOT(upSlot()));
    connect(_dn, SIGNAL(released()), this, SLOT(dnSlot()));
    connect(_saveImage, SIGNAL(released()), this, SLOT(saveImageSlot()));
    connect(_pauseButton, SIGNAL(toggled(bool)), this, SLOT(pauseSlot(bool)));
    connect(_windowButton, SIGNAL(toggled(bool)), this, SLOT(windowSlot(bool)));
    connect(_gateNumber, SIGNAL(activated(int)), this, SLOT(gateChoiceSlot(int)));
    connect(_xGrid, SIGNAL(toggled(bool)), _scopePlot, SLOT(enableXgrid(bool)));
    connect(_yGrid, SIGNAL(toggled(bool)), _scopePlot, SLOT(enableYgrid(bool)));
    connect(_blockSizeCombo, SIGNAL(activated(int)), this, SLOT(blockSizeSlot(int)));
    connect(channelButtonGroup, SIGNAL(buttonReleased(int)), this, SLOT(channelSlot(int)));

    // set the checkbox selections
    _pauseButton->setChecked(false);
    _xGrid->setChecked(true);
    _yGrid->setChecked(true);

    // initialize the book keeping for the plots.
    // This also sets up the radio buttons
    // in the plot type tab widget
    //initPlots();

    _gainKnob->setRange(-7, 7);
    _gainKnob->setTitle("Gain");

    // set the minor ticks
    _gainKnob->setScaleMaxMajor(5);
    _gainKnob->setScaleMaxMinor(5);

    _xyGraphRange = 1;
    _xyGraphCenter = 0.0;
    _knobGain = 0.0;
    _knobOffset = 0.0;
    _specGraphRange = 120.0;
    _specGraphCenter = -40.0;

    // set up the palettes
    _greenPalette = this->palette();
    _greenPalette.setColor(this->backgroundRole(), QColor("green"));
    _redPalette = _greenPalette;
    _redPalette.setColor(this->backgroundRole(), QColor("red"));

    // The initial plot type will be I and Q timeseries
    plotTypeSlot(TS_TIMESERIES_PLOT);


    // start the statistics timer
    startTimer(_statsUpdateInterval*1000);

    // let the data sources get themselves ready
    sleep(1);

}
//////////////////////////////////////////////////////////////////////
ProfilerScope::~ProfilerScope() {
}

//////////////////////////////////////////////////////////////////////
void ProfilerScope::initFFT() {

    // configure the block/fft size selection
    /// @todo add logic to insure that smallest fft size is a power of two.
    int fftSize = _config.getInt("FFTsmallest", 16);
    int maxFftSize = _config.getInt("FFTlargest", 1024);
    for (; fftSize <= maxFftSize; fftSize = fftSize*2) {
        _blockSizeChoices.push_back(fftSize);
        QString l = QString("%1").arg(fftSize);
        _blockSizeCombo->addItem(l, QVariant(fftSize));
    }
    // select the middle choice for the block size
    _blockSizeIndex = (_blockSizeChoices.size()-1)/2;
    if (_blockSizeIndex < 0)
        _blockSizeIndex = 0;
    _blockSizeCombo->setCurrentIndex(_blockSizeIndex);

    //  set up fft for power calculations:
    _fftwData.resize(_blockSizeChoices.size());
    _fftwPlan.resize(_blockSizeChoices.size());
    for (unsigned int i = 0; i < _blockSizeChoices.size(); i++) {
        // allocate the data space for fftw
        int blockSize = _blockSizeChoices[i];
        _fftwData[i] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)
                * blockSize);
        // create the plan.
        _fftwPlan[i] = fftw_plan_dft_1d(blockSize, _fftwData[i], _fftwData[i],
        FFTW_FORWARD,
        FFTW_ESTIMATE);
    }

    // create the hamming coefficients
    hammingSetup();
}

//////////////////////////////////////////////////////////////////////
void ProfilerScope::timeSeriesSlot(
        std::vector<double> I, std::vector<double> Q, double sampleRateHz,
        double tuningFreqHz) {
    if (_paused)
        return;

    processTimeSeries(I, Q);
}

//////////////////////////////////////////////////////////////////////
void ProfilerScope::saveImageSlot() {
    QString f = _config.getString("imageSaveDirectory", "c:/").c_str();

    QFileDialog d( this, tr("Save ProfilerScope Image"), f,
            tr("PNG files (*.png);;All files (*.*)"));
    d.setFileMode(QFileDialog::AnyFile);
    d.setViewMode(QFileDialog::Detail);
    d.setAcceptMode(QFileDialog::AcceptSave);
    d.setConfirmOverwrite(true);
    d.setDefaultSuffix("png");
    d.setDirectory(f);

    f = "ProfilerScope-";
    f += QDateTime::currentDateTime().toString("yyyy-MM-dd-hh-mm-ss");
    f += ".png";
    d.selectFile(f);
    if (d.exec()) {
        QStringList saveNames = d.selectedFiles();
//        _scopePlot->saveImageToFile(saveNames[0].toStdString());
        f = d.directory().absolutePath();
        _config.setString("imageSaveDirectory", f.toStdString());
    }
}
//////////////////////////////////////////////////////////////////////
void ProfilerScope::processTimeSeries(
        std::vector<double>& Idata, std::vector<double>& Qdata) {
    if (!_timeSeriesPlot)
        return;

    PlotInfo* pi = &_tsPlotInfo[_tsPlotType];
    switch (pi->getDisplayType()) {
    case ScopePlot::SPECTRUM: {

        // compute the power spectrum
        _zeroMoment = powerSpectrum(Idata, Qdata);
        displayData();
        break;
    }
    case SCOPE_PLOT_TIMESERIES:
    case SCOPE_PLOT_IVSQ: {
        I.resize(Idata.size());
        Q.resize(Qdata.size());
        I = Idata;
        Q = Qdata;
        _zeroMoment = zeroMomentFromTimeSeries(I, Q);
        displayData();
        break;
    }
    default:
        // ignore others
        break;
    }
}

//////////////////////////////////////////////////////////////////////
void ProfilerScope::displayData() {
    double yBottom = _xyGraphCenter - _xyGraphRange;
    double yTop = _xyGraphCenter + _xyGraphRange;

    QString l = QString("%1").arg(_zeroMoment, 6, 'f', 1);
    _powerDB->setText(l);

    // Time series data display
    PlotInfo* pi = &_tsPlotInfo[_tsPlotType];

    std::string xlabel;
    ScopePlot::PLOTTYPE displayType =
            (ScopePlot::PLOTTYPE) pi->getDisplayType();
    switch (displayType) {
    case ScopePlot::TIMESERIES:
        if (pi->autoscale()) {
            autoScale(I, Q, displayType);
            pi->autoscale(false);
        }
        xlabel = std::string("Time");
//        _scopePlot->TimeSeries(I, Q, yBottom, yTop, 1, xlabel, "I - Q");
        break;
    case ScopePlot::IVSQ:
        if (pi->autoscale()) {
            autoScale(I, Q, displayType);
            pi->autoscale(false);
        }
//        _scopePlot->IvsQ(I, Q, yBottom, yTop, 1, "I", "Q");
        break;
    case ScopePlot::SPECTRUM:
        if (pi->autoscale()) {
            autoScale(_spectrum, displayType);
            pi->autoscale(false);
        }
//        _scopePlot->Spectrum(_spectrum, _specGraphCenter-_specGraphRange
//                /2.0, _specGraphCenter+_specGraphRange/2.0, 1000000, false,
//                "Frequency (Hz)", "Power (dB)");
        break;
    case ScopePlot::PRODUCT:
        // include just to quiet compiler warnings
        break;
    }
}

//////////////////////////////////////////////////////////////////////
double ProfilerScope::powerSpectrum(
        std::vector<double>& Idata, std::vector<double>& Qdata) {

    int blockSize = _blockSizeChoices[_blockSizeIndex];

    _spectrum.resize(blockSize);
    int n = Idata.size();
    if (blockSize < n) {
        n = blockSize;
    }
    for (int j = 0; j < n; j++) {
        // transfer the data to the fftw input space
        _fftwData[_blockSizeIndex][j][0] = Idata[j];
        _fftwData[_blockSizeIndex][j][1] = Qdata[j];
    }
    // zero pad, if we are looking at along beam data.
    for (int j = n; j < blockSize; j++) {
        _fftwData[_blockSizeIndex][j][0] = 0;
        _fftwData[_blockSizeIndex][j][1] = 0;
    }

    // apply the hamming window to the time series
    if (_doHamming)
        doHamming();

    // caclulate the fft
    fftw_execute(_fftwPlan[_blockSizeIndex]);

    double zeroMoment = 0.0;

    // reorder and copy the results into _spectrum
    for (int i = 0; i < blockSize/2; i++) {
        double pow = _fftwData[_blockSizeIndex][i][0]
                * _fftwData[_blockSizeIndex][i][0]
                + _fftwData[_blockSizeIndex][i][1]
                        * _fftwData[_blockSizeIndex][i][1];

        zeroMoment += pow;

        pow /= blockSize*blockSize;
        pow = 10.0*log10(pow);
        _spectrum[i+blockSize/2] = pow;
    }

    for (int i = blockSize/2; i < blockSize; i++) {
        double pow = _fftwData[_blockSizeIndex][i][0]
                * _fftwData[_blockSizeIndex][i][0]
                + _fftwData[_blockSizeIndex][i][1]
                        * _fftwData[_blockSizeIndex][i][1];

        zeroMoment += pow;

        pow /= blockSize*blockSize;
        pow = 10.0*log10(pow);
        _spectrum[i - blockSize/2] = pow;
    }

    zeroMoment /= blockSize*blockSize;
    zeroMoment = 10.0*log10(zeroMoment);

    return zeroMoment;
}

////////////////////////////////////////////////////////////////////
void ProfilerScope::plotTypeSlot(
        int plotType) {

    // find out the index of the current page
    int pageNum = _typeTab->currentIndex();

    // get the radio button id of the currently selected button
    // on that page.
    int ptype = _tabButtonGroups[pageNum]->checkedId();

    // change to a raw plot type
    TS_PLOT_TYPES tstype = (TS_PLOT_TYPES)ptype;
    plotTypeChange( &_tsPlotInfo[tstype], tstype);
}

//////////////////////////////////////////////////////////////////////
void ProfilerScope::tabChangeSlot(
        QWidget* w) {
    // find out the index of the current page
    int pageNum = _typeTab->currentIndex();

    // get the radio button id of the currently selected button
    // on that page.
    int ptype = _tabButtonGroups[pageNum]->checkedId();

    // change to a raw plot type
    TS_PLOT_TYPES plotType = (TS_PLOT_TYPES)ptype;
    plotTypeChange( &_tsPlotInfo[plotType], plotType);
}

////////////////////////////////////////////////////////////////////
void ProfilerScope::plotTypeChange(
        PlotInfo* pi, TS_PLOT_TYPES newPlotType) {

    // save the gain and offset of the current plot type
    PlotInfo* currentPi;
    currentPi = &_tsPlotInfo[_tsPlotType];
    currentPi->setGain(pi->getGainMin(), pi->getGainMax(), _knobGain);
    currentPi->setOffset(pi->getOffsetMin(), pi->getOffsetMax(), _xyGraphCenter);

    // restore gain and offset for new plot type
    gainChangeSlot(pi->getGainCurrent());
    _xyGraphCenter = pi->getOffsetCurrent();

    // set the knobs for the new plot type
    _gainKnob->setValue(_knobGain);

     _tsPlotType = newPlotType;

}

////////////////////////////////////////////////////////////////////
void ProfilerScope::initPlots() {

    _pulsePlots.insert(TS_TIMESERIES_PLOT);
    _pulsePlots.insert(TS_IVSQ_PLOT);
    _pulsePlots.insert(TS_SPECTRUM_PLOT);

    _tsPlotInfo[TS_TIMESERIES_PLOT] = PlotInfo(TS_TIMESERIES_PLOT,
            SCOPE_PLOT_TIMESERIES, "I and Q", "S:  I and Q", -5.0, 5.0, 0.0,
            -5.0, 5.0, 0.0);
    _tsPlotInfo[TS_IVSQ_PLOT] = PlotInfo(TS_IVSQ_PLOT, SCOPE_PLOT_IVSQ,
            "I vs Q", "S:  I vs Q", -5.0, 5.0, 0.0, -5.0, 5.0, 0.0);
    _tsPlotInfo[TS_SPECTRUM_PLOT] = PlotInfo(TS_SPECTRUM_PLOT,
            SCOPE_PLOT_SPECTRUM, "Power Spectrum", "S:  Power Spectrum", -5.0,
            5.0, 0.0, -5.0, 5.0, 0.0);

    // remove the one tab that was put there by designer
    _typeTab->removeTab(0);

    // add tabs, and save the button group for
    // for each tab.
    QButtonGroup* pGroup;

    pGroup = addTSTypeTab("I & Q", _pulsePlots);
    _tabButtonGroups.push_back(pGroup);

    connect(_typeTab, SIGNAL(currentChanged(QWidget *)),
            this, SLOT(tabChangeSlot(QWidget*)));
}

//////////////////////////////////////////////////////////////////////
QButtonGroup* ProfilerScope::addTSTypeTab(
        std::string tabName, std::set<TS_PLOT_TYPES> types) {
    // The page that will be added to the tab widget
    QWidget* pPage = new QWidget;
    // the layout manager for the page, will contain the buttons
    QVBoxLayout* pVbox = new QVBoxLayout;
    // the button group manager, which has nothing to do with rendering
    QButtonGroup* pGroup = new QButtonGroup;

    std::set<TS_PLOT_TYPES>::iterator i;

    for (i = types.begin(); i != types.end(); i++) {
        // create the radio button
        int id = _tsPlotInfo[*i].getId();
        QRadioButton* pRadio = new QRadioButton;
        const QString label = _tsPlotInfo[*i].getLongName().c_str();
        pRadio->setText(label);

        // put the button in the button group
        pGroup->addButton(pRadio, id);
        // assign the button to the layout manager
        pVbox->addWidget(pRadio);

        // set the first radio button of the group
        // to be selected.
        if (i == types.begin()) {
            pRadio->setChecked(true);
        }
    }
    // associate the layout manager with the page
    pPage->setLayout(pVbox);

    // put the page on the tab
    _typeTab->insertTab(-1, pPage, tabName.c_str());

    // connect the button released signal to our plot type change slot.
    connect(pGroup, SIGNAL(buttonReleased(int)), this, SLOT(plotTypeSlot(int)));

    return pGroup;
}

//////////////////////////////////////////////////////////////////////
void ProfilerScope::timerEvent(
        QTimerEvent*) {

    //int rate[3];
    //for (int i = 0; i < 3; i++) {
    //		rate[i] = (_pulseCount[i] - _prevPulseCount[i])/(double)_statsUpdateInterval;
    //		_prevPulseCount[i] = _pulseCount[i];
    //}
    //	_chan0pulseCount->setNum(_pulseCount[0]/1000);
    //	_chan0pulseRate->setNum(rate[0]);
    //	_chan0errors->setNum(_errorCount[0]);
    //	_chan1pulseCount->setNum(_pulseCount[1]/1000);
    //	_chan1pulseRate->setNum(rate[1]);
    //	_chan1errors->setNum(_errorCount[1]);
    //	_chan2pulseCount->setNum(_pulseCount[2]/1000);
    //	_chan2pulseRate->setNum(rate[2]);
    //	_chan2errors->setNum(_errorCount[2]);

    //	std::cout << "Packet errors = " <<
    //	  _errorCount[0] << " " << _errorCount[1] << " " <<
    //	  _errorCount[2] << std::endl;

}

//////////////////////////////////////////////////////////////////////
void ProfilerScope::gainChangeSlot(
        double gain) {

    // keep a local copy of the gain knob value
    _knobGain = gain;

    _specGraphRange = pow(10.0, gain+2.0);

    _xyGraphRange = pow(10.0, -gain);

    _gainKnob->setValue(gain);

}

//////////////////////////////////////////////////////////////////////
void ProfilerScope::upSlot() {
    bool spectrum = false;

    if (_timeSeriesPlot) {
        PlotInfo* pi = &_tsPlotInfo[_tsPlotType];
        if (pi->getDisplayType() == ScopePlot::SPECTRUM) {
            spectrum = true;
        }
    }

    if (!spectrum) {
        _xyGraphCenter -= 0.03*_xyGraphRange;
    } else {
        _specGraphCenter -= 0.03*_specGraphRange;
    }
    displayData();
}

//////////////////////////////////////////////////////////////////////
void ProfilerScope::dnSlot() {

    bool spectrum = false;

    if (_timeSeriesPlot) {
        PlotInfo* pi = &_tsPlotInfo[_tsPlotType];
        if (pi->getDisplayType() == ScopePlot::SPECTRUM) {
            spectrum = true;
        }
    }

    if (!spectrum) {
        _xyGraphCenter += 0.03*_xyGraphRange;
    } else {
        _specGraphCenter += 0.03*_specGraphRange;
    }

    displayData();
}

//////////////////////////////////////////////////////////////////////
void ProfilerScope::autoScale(
        std::vector<double>& data, ScopePlot::PLOTTYPE displayType) {
    if (data.size() == 0)
        return;

    // find the min and max
    double min = *std::min_element(data.begin(), data.end());
    double max = *std::max_element(data.begin(), data.end());

    // adjust the gains
    adjustGainOffset(min, max, displayType);
}

//////////////////////////////////////////////////////////////////////
void ProfilerScope::autoScale(
        std::vector<double>& data1, std::vector<double>& data2,
        ScopePlot::PLOTTYPE displayType) {
    if (data1.size() == 0 || data2.size() == 0)
        return;

    // find the min and max
    double min1 = *std::min_element(data1.begin(), data1.end());
    double min2 = *std::min_element(data2.begin(), data2.end());
    double min = std::min(min1, min2);

    double max1 = *std::max_element(data1.begin(), data1.end());
    double max2 = *std::max_element(data2.begin(), data2.end());
    double max = std::max(max1, max2);

    // adjust the gains
    adjustGainOffset(min, max, displayType);

}

//////////////////////////////////////////////////////////////////////
void ProfilerScope::adjustGainOffset(
        double min, double max, ScopePlot::PLOTTYPE displayType) {
    if (displayType == ScopePlot::SPECTRUM) {
        // currently in spectrum plot mode
        _specGraphCenter = min + (max-min)/2.0;
        _specGraphRange = 3*(max-min);
        _knobGain = -log10(_specGraphRange);
    } else {
        double factor = 0.8;
        _xyGraphCenter = (min+max)/2.0;
        _xyGraphRange = (1/factor)*(max - min)/2.0;
        if (min == max ||
        		isnan(min) ||
        		isnan(max) ||
        		isinf(min) ||
        		isinf(max))
        	_xyGraphRange = 1.0;
        //std::cout << "min:"<<min<<"  max:"<<max<<"     _xxGraphRange is " << _xyGraphRange << "\n";
        _knobGain = -log10(_xyGraphRange);
        _gainKnob->setValue(_knobGain);
    }
}


//////////////////////////////////////////////////////////////////////
void
ProfilerScope::newTSItemSlot(ProfilerDDS::TimeSeries* pItem) {
	emit returnTSItem(pItem);
}
//////////////////////////////////////////////////////////////////////
void ProfilerScope::autoScaleSlot() {
    PlotInfo* pi;

    pi = &_tsPlotInfo[_tsPlotType];

    pi->autoscale(true);
}

//////////////////////////////////////////////////////////////////////
void ProfilerScope::pauseSlot(
        bool p) {
    _paused = p;
}

//////////////////////////////////////////////////////////////////////
void ProfilerScope::channelSlot(
        int c) {
    _channel = c;
}

//////////////////////////////////////////////////////////////////////
void ProfilerScope::gateChoiceSlot(
        int index) {
    _gateChoice = _gates[index];

}

//////////////////////////////////////////////////////////////////////
void ProfilerScope::blockSizeSlot(
        int index) {

	_blockSizeIndex = index;

    // recalculate the hamming coefficients. _blockSizeIndex
	// must be set correctly before calling this
    hammingSetup();

}

////////////////////////////////////////////////////////////////////////

double ProfilerScope::zeroMomentFromTimeSeries(
        std::vector<double>& I, std::vector<double>& Q) {
    double p = 0;
    int n = I.size();

    for (unsigned int i = 0; i < I.size(); i++) {
        p += I[i]*I[i] + Q[i]*Q[i];
    }

    p /= n;
    p = 10.0*log10(p);
    return p;
}

////////////////////////////////////////////////////////////////////////
void
ProfilerScope::doHamming() {

  int blockSize = _blockSizeChoices[_blockSizeIndex];

  for (int i = 0; i < blockSize; i++) {
    _fftwData[_blockSizeIndex][i][0] *= _hammingCoefs[i];
    _fftwData[_blockSizeIndex][i][1] *= _hammingCoefs[i];
  }
}
////////////////////////////////////////////////////////////////////////

void
ProfilerScope::hammingSetup() {

   int blockSize = _blockSizeChoices[_blockSizeIndex];

  _hammingCoefs.resize(blockSize);

  for (int i = 0; i < blockSize; i++) {
    _hammingCoefs[i] = 0.54 - 0.46*(cos(2.0*M_PI*i/(blockSize-1)));
  }

}

////////////////////////////////////////////////////////////////////////

void
ProfilerScope::windowSlot(bool flag) {
	_doHamming = flag;
}
