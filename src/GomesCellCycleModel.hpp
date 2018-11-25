#ifndef GOMESCELLCYCLEMODEL_HPP_
#define GOMESCELLCYCLEMODEL_HPP_

#include "AbstractSimpleCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"
#include "Cell.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "GomesRetinalNeuralFates.hpp"
#include "SmartPointers.hpp"
#include "ColumnDataWriter.hpp"
#include "LogFile.hpp"
#include "CellLabel.hpp"

/*******************************************
 * GOMES CELL CYCLE MODEL
 * As described in Gomes et al. 2011 [Gomes2011] doi: 10.1242/dev.059683
 *
 * USE: By default, GomesCellCycleModels are constructed with the parameter fit reported in [Gomes2011].
 * Cell cycle length, mitotic mode, and postmitotic fate of cells are determined by independent random variables
 * PP = symmetric proliferative mitotic mode, both progeny remain mitotic
 * PD = asymmetric proliferative mitotic mode, one progeny exits the cell cycle and differentiates
 * DD = symmetric differentiative mitotic mode, both progeny exit the cell cycle and differentiate
 *
 * Change default model parameters with SetModelParameters(<params>);
 * Set AbstractCellProperties for differentiated neural types with SetModelProperties();
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
 **********************************************/

class GomesCellCycleModel : public AbstractSimpleCellCycleModel
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
    void WriteDebugData(double percentile);

protected:
    //mode/output variables
    bool mOutput;
    double mEventStartTime;
    bool mSequenceSampler;
    bool mSeqSamplerLabelSister;
    //debug writer stuff
    bool mDebug;
    int mTimeID;
    std::vector<int> mVarIDs;
    boost::shared_ptr<ColumnDataWriter> mDebugWriter;
    //model parameters and state memory vars
    double mNormalMu;
    double mNormalSigma;
    double mPP;
    double mPD;
    double mpBC;
    double mpAC;
    double mpMG;
    unsigned mMitoticMode;
    unsigned mSeed;
    boost::shared_ptr<AbstractCellProperty> mp_PostMitoticType;
    boost::shared_ptr<AbstractCellProperty> mp_RPh_Type;
    boost::shared_ptr<AbstractCellProperty> mp_BC_Type;
    boost::shared_ptr<AbstractCellProperty> mp_AC_Type;
    boost::shared_ptr<AbstractCellProperty> mp_MG_Type;
    boost::shared_ptr<AbstractCellProperty> mp_label_Type;

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
    GomesCellCycleModel(const GomesCellCycleModel& rModel);

public:

    /**
     * Constructor - just a default, mBirthTime is set in the AbstractCellCycleModel class.
     */
    GomesCellCycleModel();

    /**
     * SetCellCycleDuration() method to set length of cell cycle (lognormal distribution as specified in [Gomes2011])
     */
    void SetCellCycleDuration();

    /**
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

    /** Overridden InitialiseDaughterCell() method. Used to implement asymmetric mitotic mode*/
    void InitialiseDaughterCell();

    /* Model setup functions
     * Set lognormal cell cycle curve properties with *mean and std of corresponding NORMAL curve*
     * Model requires valid AbstractCellProperties to assign postmitotic fate;
     * Defaults are found in GomesRetinalNeuralFates.hpp
     * (RodPohotreceptor, AmacrineCell, BipolarCell, MullerGlia)
     */

    void SetModelParameters(const double normalMu = 3.9716, const double normalSigma = .32839, const double PP = .055,
                            const double PD = .221, const double pBC = .128, const double pAC = .106, const double pMG =
                                    .028);
    void SetModelProperties(boost::shared_ptr<AbstractCellProperty> p_RPh_Type,
                            boost::shared_ptr<AbstractCellProperty> p_AC_Type,
                            boost::shared_ptr<AbstractCellProperty> p_BC_Type,
                            boost::shared_ptr<AbstractCellProperty> p_MG_Type);

    //This should normally be a DifferentiatedCellProliferativeType
    void SetPostMitoticType(boost::shared_ptr<AbstractCellProperty> p_PostMitoticType);

    //Functions to enable per-cell mitotic mode logging for mode rate & sequence sampling fixtures
    //Uses singleton logfile
    void EnableModeEventOutput(double eventStart, unsigned seed);
    void EnableSequenceSampler(boost::shared_ptr<AbstractCellProperty> label);

    //More detailed debug output. Needs a ColumnDataWriter passed to it
    //Only declare ColumnDataWriter directory, filename, etc; do not set up otherwise
    void EnableModelDebugOutput(boost::shared_ptr<ColumnDataWriter> debugWriter);

    //Not used, but must be overwritten lest GomesCellCycleModels be abstract
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
CHASTE_CLASS_EXPORT(GomesCellCycleModel)

#endif /*GOMESCELLCYCLEMODEL_HPP_*/
