#ifndef PROFILERSCOPE_H_
#define PROFILERSCOPE_H_

#include <QWidget>
#include <QPalette>
#include <QButtonGroup>

#include <qevent.h>
#include <deque>
#include <set>
#include <map>
#include <fftw3.h>

// Components from the QtToolbox
#include "ScopePlot.h"
#include "Knob.h"
#include "QtConfig.h"

// The designer generated header file.
#include "ui_AScope.h"

// PlotInfo knows the characteristics of a plot
#include "PlotInfo.h"

/**
 AScope provides a traditional real-time Ascope display of
 eldora time series data and computed products. It is implmented
 with Qt, and uses the QtToolbox::ScopePlot as the primary display.
 I&Q, I versus Q, IQ power spectrum, and computed product displays
 are user selectable. The data can be displayed either along the beam
 for all gates, or in time for one gate. Users may select the
 fft block size and the gate to be displayed.

 AScope is simply a data consumer; it does not know
 anything about the data provider. Signals and slots
 used to coordinate with other components.

 It is the responsibility of the data provider to feed data
 at a desired rate. AScope will attempt to render all data
 delivered.
 **/
class AScope : public QWidget, private Ui::AScope {
    Q_OBJECT
        /// types of plots available in the scope plot.
        enum SCOPE_PLOT_TYPES {
            SCOPE_PLOT_TIMESERIES,
            SCOPE_PLOT_IVSQ,
            SCOPE_PLOT_SPECTRUM
        };

        /// Time series plot types.
        enum TS_PLOT_TYPES {
            TS_TIMESERIES_PLOT, ///<  time series I and Q plot
            TS_IVSQ_PLOT,       ///<  time series I versus Q plot
            TS_SPECTRUM_PLOT    ///<  time series power spectrum plot
        };
        
     public:
        /// The timeseries type for importing data. The actual data
        /// are passed by reference, hopefully eliminating an
        /// unnecessary copy.
        typedef struct {
        	/// I and Q for each beam is in a vector containing I,Q for each gate.
        	/// IQbeams contains pointers to each IQ vector for all
        	/// of the beams in the timeseries. The length of the timeseries
        	/// can be found from IQbeams.size().
        	std::vector<float*> IQbeams;
        	/// The number of gates
        	int gates;
        	/// The channel id
        	int chanId;
        	/// An opaque pointer that can be used to store
        	/// anything that the caller wants to track along 
        	/// with the TimeSeries. This will be useful when 
        	/// the TimeSeries is returned to the orginator, 
        	/// if for example an associated object such as a
        	/// DDS sample needs to be returned to DDS.
        	void* handle;
        } TimeSeries;

        AScope(
                QWidget* parent = 0);
        virtual ~AScope();

    signals:
		/// emit this signal to alert the client that we
		/// are finished with this item. AScope::TimeSeries
		/// contains an opaque handle that the client can
		/// use to keep track of this item between the 
		/// triggering of newTSItemSlot() and the emmiting
		/// of returnTSItem().
		void returnTSItem(AScope::TimeSeries pItem);

    public slots:
		/// Feed new timeseries data via this slot.
		/// @param pItem This contains some metadata and pointers to I/Q data
		void newTSItemSlot(AScope::TimeSeries pItem);
       /// Call when the plot type is changed. This function
        /// must determine which of the two families of
        /// plots, _tsPlotInfo, or _productPlotInfo, the
        /// previous and new plot types belong to.
        virtual void plotTypeSlot(
                int plotType);
        /// call to save the current plotting parameters for the
        /// current plot type, and reload the parameters for the
        /// the new plot type.
        void plotTypeChange(
                PlotInfo* pi,
                    TS_PLOT_TYPES plotType);
        /// A different tab has been selected. Change the plot type to the
        /// currently selected button on that tab.
        void tabChangeSlot(
                QWidget* w);
        /// The gain knob value has changed.
        virtual void gainChangeSlot(
                double);
        /// slide the plot up.
        virtual void upSlot();
        /// Slide the plot down.
        virtual void dnSlot();
        /// Initiate an autoscale. A flag is set; during the next
        /// pulse reception an autoscale computation is made.
        virtual void autoScaleSlot();
        /// Save the scope display to a PNG file.
        void saveImageSlot();
        /// Pause the plotting. Any received data are ignored.
        /// @param p True to enable pause.
        void pauseSlot(
                bool p);
        /// Select the channel
        /// @param c The channel (0-3)
        void channelSlot(
                int c);
        /// Select the gate
        /// @param index The index from the combo box of the selected gate.
        void gateChoiceSlot(
                int index);
        /// Select the block size
        /// @param size The block size. It must be a power of two.
        void blockSizeSlot(
                int size);
        /// Enable/disable windowing
        void windowSlot(bool);

    protected:
        /// Initialize the fft calculations. The minimum
        /// size is 8. The max size is time series length
        /// to the largest power of 2.
        /// Fftw plans and data arrays are allocated for all
        /// powers of two within this range. The block size
        /// combo selector is initialized.
        void initFFT(int tsLength);
        /// Initialize the gate selection 
        /// @param gates The number of gates
        void initGates(int gates);
        /// Initialize the channel selection
        /// @param channels The number of channels
		void initChans(int channels);
        /// Emit a signal announcing the desired gate mode,
        /// either along beam, or one gate. The channel select,
        /// gate choice and (for one gate mode) data block
        /// size will be part of the emitted signal.
        void dataMode();
        /// Send the data for the current plot type to the ScopePlot.
        void displayData();
        /// setup the hamming coefficients
        void hammingSetup();
        /// apply the hamming filter
        void doHamming();
        ///	cumulative error count
        int _errorCount[3];
        ///  last pulse number
        long long _lastPulseNum[3];
        double _knobGain;
        double _knobOffset;
        double _xyGraphRange;
        double _xyGraphCenter;
        double _specGraphRange;
        double _specGraphCenter;
        /// Autoscale based on a set of data.
        /// @param data The data series to be analyzed.
        /// @param displayType The ScopePlot::PLOTTYPE of the display
        void autoScale(
                std::vector<double>& data,
                ScopePlot::PLOTTYPE displayType);
        /// Autoscale based on two sets of data.
        /// @param data1 The first data series to be analyzed.
        /// @param data2 The second data series to be analyzed.
        /// @param displayType The ScopePlot::PLOTTYPE of the display
        void autoScale(
                std::vector<double>& data1,
                    std::vector<double>& data2,
                    ScopePlot::PLOTTYPE displayType);
        /// Initialize the combo box choices and FFTs.
        /// @param channels The number of channels,
        /// @param tsLength The time series length
        /// @param gates The number of gates
        void initCombos(int channels, int tsLength, int gates);
        /// Adjust the _graphRange and _graphOffset values.
        /// @param min Desired scale minimum
        /// @param max Desired scale maximum
        /// @param displayType The ScopePlot::PLOTTYPE of the display
        void adjustGainOffset(
                double min,
                    double max,
                    ScopePlot::PLOTTYPE displayType);
        /// save the button group for each tab,
        /// so that we can find the selected button
        /// and change the plot type when tabs are switched.
        std::vector<QButtonGroup*> _tabButtonGroups;
        /// This set contains PLOTTYPEs for all raw data plots
        std::set<TS_PLOT_TYPES> _pulsePlots;
        /// Holds I data to display for time series and I vs. Q
        std::vector<double> I;
        /// Holds Q data to display for time series and I vs. Q display
        std::vector<double> Q;
        /// Holds power spectrum values for display.
        std::vector<double> _spectrum;
        // how often to update the statistics (in seconds)
        int _statsUpdateInterval;
        /// Set true if time series plots are selected, false for product type plots
        bool _timeSeriesPlot;
        /// The current selected plot type.
        TS_PLOT_TYPES _tsPlotType;
        // The builtin timer will be used to calculate beam statistics.
        void timerEvent(
                QTimerEvent*);
        /// The hamming window coefficients
        std::vector<double> _hammingCoefs;
        /// The index of the current fft/block size
        /// selection in _blockSizeChoices
        int _blockSizeIndex;
        /// The possible block/fftw size choices.
        std::vector<int> _blockSizeChoices;
        ///	The fftw plan. This is a handle used by
        ///	the fftw routines.
        std::vector<fftw_plan> _fftwPlan;
        ///	The fftw data array. The fft will
        //	be performed in place, so both input data
        ///	and results are stored here.
        std::vector<fftw_complex*> _fftwData;
        //	power correction factor applied to (uncorrected) powerSpectrum() output
        double _powerCorrection;
        /// Set true if the Hamming window should be applied
        bool _doHamming;
        /// Process time series data.
        /// @param Idata The I values
        /// @param Qdata The Q values
        void processTimeSeries(
                std::vector<double>& Idata,
                    std::vector<double>& Qdata);
        /// Compute the power spectrum. The input values will come
        /// I[]and Q[], the power spectrum will be written to
        /// _spectrum[]
        /// @param Idata The I time series.
        /// @param Qdata The Q time series.
        /// @return The zero moment
        double powerSpectrum(
                std::vector<double>& Idata,
                    std::vector<double>& Qdata);
        /// For each TS_PLOT_TYPES, there will be an entry in this map.
        std::map<TS_PLOT_TYPES, PlotInfo> _tsPlotInfo;
        /// This set contains PLOTTYPEs for all timeseries plots
        std::set<TS_PLOT_TYPES> _timeSeriesPlots;
        /// initialize all of the book keeping structures
        /// for the various plots.
        void initPlots();
        /// add a ts tab to the plot type selection tab widget.
        /// Radio buttons are created for all of specified
        /// plty types, and grouped into one button group.
        /// _tsPlotInfo provides the label information for
        /// the radio buttons.
        /// @param tabName The title for the tab.
        /// @param types A set of the desired TS_PLOT_TYPES types
        /// @return The button group that the inserted buttons
        /// belong to.
        QButtonGroup* addTSTypeTab(
                std::string tabName,
                    std::set<TS_PLOT_TYPES> types);
        /// Calculate the zeroth moment, using the time
        /// series for input.
        double zeroMomentFromTimeSeries(
                std::vector<double>& I,
                    std::vector<double>& Q);
        /// The configuration for AScope
        QtConfig _config;
        /// The button group for channel selection
        QButtonGroup* _chanButtonGroup;
        /// Palette for making the leds green
        QPalette _greenPalette;
        /// Platette for making the leds red
        QPalette _redPalette;
        /// Set true if the plot graphics are paused
        bool _paused;
        /// The signal power, computed directly from the I&Q
        /// data, or from the power spectrum
        double _zeroMoment;
        /// The choice of channels (0-3)
        int _channel;
        /// The selected gate, zero based.
        int _gateChoice;
        /// Set false to cause initialization of blocksize and 
        /// gate choices when the first data is received.
        bool _combosInitialized;
       
};


#endif /*PROFSCOPE_H_*/
