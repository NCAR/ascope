#ifndef PROFSCOPE_H_
#define PROFSCOPE_H_

#include <QDialog>
#include <QPalette>
#include <qevent.h>
#include <deque>
#include <set>
#include <map>
#include <fftw3.h>

// Coponents from the QtToolbox
#include "ScopePlot.h"
#include "Knob.h"
#include "QtConfig.h"

// The designer generated header file.
#include "ui_EldoraScope.h"

// PlotInfo knows the characteristics of a plot
#include "PlotInfo.h"

// Symbolic names for product types
#include "ProductTypes.h"

/** 
 EldoraScope provides a traditional real-time Ascope display of 
 eldora time series data and computed products. It is implmented
 with Qt, and uses the QtToolbox::ScopePlot as the primary display.
 I&Q, I versus Q, IQ power spectrum, and computed product displays
 are user selectable. The data can be displayed either along the beam
 for all gates, or in time for one gate. Users may select the 
 fft block size and the gate to be displayed.
 
 EldoraScope is simply a data consumer; it does not know
 anything about the data provider. Signals and slots
 used to coordinate with other components. Data are expected
 to be delivered in a mode matching the current display
 mode. EldoraScope announces the mode by emitting either
 an alongBeam signal or a oneGate signal. Data are then delivered
 to EldoraScope by calling the newTimeSeriesSlot() and newProductSlot().
 
 It is the responsibility of the data provider to feed data
 at a desired rate. EldoraScope will attempt to render all data
 delivered to newTimeSeriesSlot() and newProductSlot().

 EldoraScope is configured via EldoraScope.ini.
 **/
class EldoraScope : public QDialog, public Ui::EldoraScope {
    Q_OBJECT
        /// The display can either show all gates along a beam, or
        /// values in time for a selected gate.
        enum GATE_MODE {
            ALONG_BEAM, ///< Display all gates along a beam
            ONE_GATE ///< Display values in time for a selected gate
        };

        /// types of plots available in the scope plot.
        enum SCOPE_PLOT_TYPES {
            SCOPE_PLOT_TIMESERIES,
            SCOPE_PLOT_IVSQ,
            SCOPE_PLOT_SPECTRUM,
            SCOPE_PLOT_PRODUCT
        };

        /// Time series plot types.
        enum TS_PLOT_TYPES {
            TS_TIMESERIES_PLOT, ///<  time series I and Q plot
            TS_IVSQ_PLOT, ///<  time series I versus Q plot
            TS_SPECTRUM_PLOT ///<  time series power spectrum plot 
        };

     public:
        EldoraScope(
                QDialog* parent = 0);
        virtual ~EldoraScope();

    signals:

    /// Emitted to announce that TS data should be delivered 
    /// in ONE_GATE mode.
    /// @param channel The selected channel
    /// @param gate The selected gate
    /// @param n The number of points to deliver for the selected gate
    void oneGateTSSignal(
            int channel,
                bool forwardRadar,
                int gate,
                int n);
    /// emmited to indicate that TS data should be delivered for 
    /// all gates along a beam
    /// @param channel The selected channel
    void alongBeamTSSignal(
            int channel,
            bool forwardRadar);
    
    /// Emitted to announce that Product data should be delivered 
    /// in ONE_GATE mode.
    /// @param product The selected product
    /// @param gate The selected gate
    /// @param n The number of points to deliver for the selected gate
    void oneGateProductSignal(
            PRODUCT_TYPES product,
                bool forwardRadar,
                int gate,
                int n);
    /// Emmited to indicate that Product data should be delivered for 
    /// all gates along a beam
    /// @param channel The selected channel
    /// @param forwardRadar Set true if the forward radar, false otherwise.
    void alongBeamProductSignal(
            PRODUCT_TYPES product,
            bool forwardRadar);
    
    public slots:
        /// Feed new timeseries data via this slot. The data 
        /// vectors must be of the same length and non-zero; otherwise they
        /// will be ignored. The vector lengths can change between calls,
        /// and the plot will respond appropriately. If the plot is currently 
        /// configured for a time series display (I&Q, IvsQ or spectrum), the
        /// new data will be displayed.
        /// @param I A vector I values
        /// @param Q A vector of Q values
        /// @param sampleRateHz The sample rate of the I/Q data, in hz.
        /// @param tuningFrequencyHz The current frequency of the receiver.
        void timeSeriesSlot(
                std::vector<double> I,
                    std::vector<double> Q,
                    double sampleRateHz,
                    double tuningFreqHz);
        /// Call when data is available on the product data socket.
        /// @param p the product data
        /// @param radarId either EldoraDDS::Fore or EldoraDDS::Aft
        /// @param elDegrees The antenna pointing elevation, degrees
        /// @param prodType The product type, from PRODUCT_TYPES
        /// @param gateSpacingMeters The width in meters of each gate
        /// @param dwellWidth The angular width of one dwell
        /// @param airSpdCorr The airspeed correction to be added to the 
        /// radial velocity, if desired.
        /// @param rollAngle The aircraft roll angle.
        /// @param nyquistVelocity The radar nyquist velocity, m/s.        
        /// @param altitudeMSL Altitude, meters abocw sea level
        /// @param latitude Latitude, degrees
        /// @param longitude Longitude, degrees

        void productSlot(std::vector<double> p, 
                         int radarId, 
                         float elDegrees, 
                         int prodType, 
                         float gateSpacingMeters,
                         double dwellWidth,
                         double airSpdCorr,
                         double rollAngle,
                         double nyquistVelocity,
                         double altitudeMSL,
                         double latitude,
                         double longitude);
        /// Call to set the list of available gates in the timeseries.
        /// @param gates A list of possible gates in the timeseries. It will be zero based. 
        /// Add firstgate to stablish the true gate number.
        /// @param firstgate The gate number of the first gate.
        void tsGateListSlot(
                std::vector<int> gates, unsigned short firstgate);
        /// Call when the plot type is changed. This function 
        /// must determine which of the two families of
        /// plots, _tsPlotInfo, or _productPlotInfo, the
        /// previous and new plot types belong to.
        virtual void plotTypeSlot(
                int plotType);
        /// call to save the current plotting parameters for the
        /// current plot type, and reload the parameters for the 
        /// the new plot type. It handles both pulse and beam
        /// displays. pulsePlot is used to differentiate between the
        /// two.
        void plotTypeChange(
                PlotInfo* pi,
                    TS_PLOT_TYPES plotType,
                    PRODUCT_TYPES prodType,
                    bool pulsePlot);
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
        /// Set the gate display mode
        /// @param m The gate mode, either ALONG_BEAM or ONE_GATE
        void gateModeSlot(
                int m);
        /// Select the channel
        /// @param c The channel (1-4)
        void channelSlot(
                int c);
        /// Select the radar
        /// @param forwardRadar True if forward radar, false if aft
        void radarSlot(
                int forwardRadar);
        /// Select the gate
        /// @param g The index from the combo box of the selected gate.
        void gateChoiceSlot(
                int index);
        /// Select the block size
        /// @param The block size. It must be a power of two.
        void blockSizeSlot(
                int);
        /// Enable/disable windowing
        void windowSlot(bool);

    protected:
        /// Initialize the fft calculations. The minimum and
        /// maximum fft sizes are read from the configuration.
        /// Fftw plans and data arrays are allocated for all
        /// powers of two within this range. The block size
        /// combo selector is initialized.
        void initFFT();
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
        /// set the _graphRange and _graphOffset based
        /// on the single data series.
        /// @param data The data series to be analyzed.
        void autoScale(
                std::vector<double>& data,
                ScopePlot::PLOTTYPE displayType);
        /// set the _graphRange and _graphOffset based
        /// on the two data series.
        /// @param data1 The first data series to be analyzed.
        /// @param data2 The second data series to be analyzed.
        void autoScale(
                std::vector<double>& data1,
                    std::vector<double>& data2,
                    ScopePlot::PLOTTYPE displayType);
        /// Adjust the _graphRange and _graphOffset values.
        /// @param min Desired scale minimum
        /// @param max Desired scale maximum
        void adjustGainOffset(
                double min,
                    double max,
                    ScopePlot::PLOTTYPE displayType);
        /// Holds I data to display for time series and I vs. Q 	
        std::vector<double> I;
        /// Holds Q data to display for time series and I vs. Q display
        std::vector<double> Q;
        /// Holds power spectrum values for display.
        std::vector<double> _spectrum;
        /// Used to collect product data from beams
        std::vector<double> _ProductData;
        // how often to update the statistics (in seconds)
        int _statsUpdateInterval;
        /// Set true if time series plots are selected, false for product type plots
        bool _timeSeriesPlot;
        /// The current selected plot type.
        TS_PLOT_TYPES _tsPlotType;
        /// The current selected product type.
        PRODUCT_TYPES _productPlotType;
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
        /// Process product data
        void processProduct(std::vector<double>& p);
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
        /// For each PRODUCT_PLOT_TYPES, there will be an entry in this map.
        std::map<PRODUCT_TYPES, PlotInfo> _productPlotInfo;
        /// This set contains PLOTTYPEs for all timeseries plots
        std::set<TS_PLOT_TYPES> _timeSeriesPlots;
        /// This set contains PLOTTYPEs for all raw data plots
        std::set<TS_PLOT_TYPES> _pulsePlots;
        /// This set contains PLOTTYPEs for all S band moments plots
        std::set<PRODUCT_TYPES> _productPlots;
        /// This set contains PLOTTYPEs for all X band moments plots
        std::set<PRODUCT_TYPES> _xMomentsPlots;
        /// save the button group for each tab,
        /// so that we can find the selected button
        /// and change the plot type when tabs are switched.
        std::vector<QButtonGroup*> _tabButtonGroups;
        /// initialize all of the book keeping structures
        /// for the various plots.
        void initPlots();
        /// add a rw plot tab to the plot type selection tab widget.
        /// Radio buttons are created for all of specified
        /// plty types, and grouped into one button group.
        /// _tsPlotInfo provides the label information for
        /// the radio buttons.
        /// @param tabName The title for the tab.
        /// @param types A set of the desired TS_PLOT_TYPES types 
        /// @return The button group that the inserted buttons
        /// belong to.
        QButtonGroup* addProductTypeTab(
                std::string tabName,
                    std::set<PRODUCT_TYPES> types);
        /// add a products tab to the plot type selection tab widget.
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
        /// The configuration for EldoraScope
        QtConfig _config;
        /// Palette for making the leds green
        QPalette _greenPalette;
        /// Platette for making the leds red
        QPalette _redPalette;
        /// Set true if the plot graphics are paused
        bool _paused;
        /// The selected gate mode
        GATE_MODE _gateMode;
        /// A list of available gates. This is set via the 
        /// gateList slot.
        std::vector<int> _gates;
        /// The choice of channels (1-4)
        int _channel;
        /// The selected gate, zero based.
        int _gateChoice;
        /// The signal power, computed directly from the I&Q
        /// data, or from the power spectrum
        double _zeroMoment;
        /// True if forward radar, false otherwise
        bool _forwardRadar;
};


#endif /*PROFSCOPE_H_*/
