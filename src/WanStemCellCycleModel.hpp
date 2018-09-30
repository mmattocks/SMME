#ifndef WANSTEMCELLCYCLEMODEL_HPP_
#define WANSTEMCELLCYCLEMODEL_HPP_

#include "AbstractSimpleCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"
#include "Cell.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "SmartPointers.hpp"
#include "ColumnDataWriter.hpp"
#include "LogFile.hpp"

#include "HeCellCycleModel.hpp"

/***********************************
 * WAN STEM CELL CYCLE MODEL
 * A simple model standing in for the CMZ "stem cells" described in'
 * Wan et al. 2016 [Wan2016]
 *
 * USE: By default, WanStemCellCycleModels consistently divide asymmetrically.
 * InitialiseDaughterCell() marks offspring for RPC fate
 * So-marked cells are given HeCellCycleModels for their next division.
 *
 * Change default model parameters with SetModelParameters(<params>);
 *
 * 2 per-model-event output modes:
 * EnableModeEventOutput() enables mitotic mode event logging-all cells will write to the singleton log file
 * EnableModelDebugOutput() enables more detailed debug output, each seed will have its own file written to
 * by a ColumnDataWriter passed to it from the test
 * (eg. by the SetupDebugOutput helper function in the project simulator)
 *
 * 1 mitotic-event-sequence sampler (only samples one "path" through the lineage):
 * EnableSequenceSampler() - one "sequence" of progenitors writes mitotic event type to a string in the singleton log file
 *
 ************************************/

class WanStemCellCycleModel : public AbstractSimpleCellCycleModel
{
    friend class TestSimpleCellCycleModels;

private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the cell-cycle model and random number generator, never used directly - boost uses this.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractSimpleCellCycleModel>(*this);

        SerializableSingleton<RandomNumberGenerator>* p_wrapper =
                RandomNumberGenerator::Instance()->GetSerializationWrapper();
        archive & p_wrapper;
        archive & mCellCycleDuration;
    }

    //Private write functions for models
    void WriteModeEventOutput();
    void WriteDebugData();

protected:
    //mode/output variables
    bool mOutput;
    double mEventStartTime;
    //debug writer stuff
    bool mDebug;
    int mTimeID;
    std::vector<int> mVarIDs;
    boost::shared_ptr<ColumnDataWriter> mDebugWriter;
    //model parameters and state memory vars
    bool mRPCfate;
    double mGammaShift;
    double mGammaShape;
    double mGammaScale;
    unsigned mSeed;
    bool mTimeDependentCycleDuration;
    double mPeakRateTime;
    double mIncreasingRateSlope;
    double mDecreasingRateSlope;
    double mBaseGammaScale;
    boost::shared_ptr<AbstractCellProperty> mp_TransitType;
    std::vector<double> mHeParamVector;

    /**
     * Protected copy-constructor for use by CreateCellCycleModel().
     *
     * The only way for external code to create a copy of a cell cycle model
     * is by calling that method, to ensure that a model of the correct subclass is created.
     * This copy-constructor helps subclasses to ensure that all member variables are correctly copied when this happens.
     *
     * This method is called by child classes to set member variables for a daughter cell upon cell division.
     * Note that the parent cell cycle model will have had ResetForDivision() called just before CreateCellCycleModel() is called,
     * so performing an exact copy of the parent is suitable behaviour. Any daughter-cell-specific initialisation
     * can be done in InitialiseDaughterCell().
     *
     * @param rModel the cell cycle model to copy.
     */
    WanStemCellCycleModel(const WanStemCellCycleModel& rModel);

public:

    /**
     * Constructor - just a default, mBirthTime is set in the AbstractCellCycleModel class.
     */
    WanStemCellCycleModel();

    /**
     * SetCellCycleDuration() method to set length of cell cycle
     */
    void SetCellCycleDuration();

    /**pro
     * Overridden builder method to create new copies of
     * this cell-cycle model.
     *
     * @return new cell-cycle model
     */
    AbstractCellCycleModel* CreateCellCycleModel();

    /**
     * Overridden ResetForDivision() method.
     * Contains general mitotic mode logic
     **/
    void ResetForDivision();

    /**
     * Overridden InitialiseDaughterCell() method.
     * Used to apply sister-cell time shifting (cell cycle duration, deterministic phase boundaries)
     * Used to implement asymmetric mitotic mode
     * */
    void InitialiseDaughterCell();

    /*Model setup functions for standard He (SetModelParameters) and deterministic alternative (SetDeterministicMode) models
     * Default parameters are from refits of He et al + deterministic alternatives
     * He 2012 params: mitoticModePhase2 = 8, mitoticModePhase3 = 15, p1PP = 1, p1PD = 0, p2PP = .2, p2PD = .4, p3PP = .2, p3PD = 0
     * gammaShift = 4, gammaShape = 2, gammaScale = 1, sisterShift = 1
     */
    void SetModelParameters(double gammaShift = 4, double gammaShape = 2, double gammaScale = 1, std::vector<double> heParamVector = { 8, 15, 1, 0, .2, .4, .2, 0, 4, 2, 1, 1 });
    void SetTimeDependentCycleDuration(double peakRateTime, double increasingSlope, double decreasingSlope);

    //This should normally be a DifferentiatedCellProliferativeType
    void SetTransitType(boost::shared_ptr<AbstractCellProperty> p_PostMitoticType);

    //Functions to enable per-cell mitotic mode logging for mode rate & sequence sampling fixtures
    //Uses singleton logfile
    void EnableModeEventOutput(double eventStart, unsigned seed);

    //More detailed debug output. Needs a ColumnDataWriter passed to it
    //Only declare ColumnDataWriter directory, filename, etc; do not set up otherwise
    void EnableModelDebugOutput(boost::shared_ptr<ColumnDataWriter> debugWriter);

    //Not used, but must be overwritten lest WanStemCellCycleModels be abstract
    double GetAverageTransitCellCycleTime();
    double GetAverageStemCellCycleTime();

    /**
     * Overridden OutputCellCycleModelParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputCellCycleModelParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(WanStemCellCycleModel)

#endif /*WANSTEMCELLCYCLEMODEL_HPP_*/
