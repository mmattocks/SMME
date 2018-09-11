#ifndef GOMESCELLCYCLEMODEL_HPP_
#define GOMESCELLCYCLEMODEL_HPP_

#include "AbstractSimpleCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "ColumnDataWriter.hpp"
#include "Cell.hpp"
#include <RandomNumberGenerator.hpp>
#include <SimulationTime.hpp>
#include "HeAth5Mo.hpp"
#include "SmartPointers.hpp"
#include "LogFile.hpp"

/**
 * A stochastic cell-cycle model where cells divide with a stochastic cell cycle duration
 * with the length of the cell cycle drawn from a uniform distribution
 * on [mMinCellCycleDuration, mMaxCellCycleDuration].
 *
 * If the cell is differentiated, then the cell cycle duration is set to be infinite,
 * so that the cell will never divide.
 */
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

protected:
    bool mDebug;
    bool mOutput;
    bool mSequenceSampler;
    bool mSeqSamplerLabelSister;
    int mTimeID;
    std::vector<int> mVarIDs;
    boost::shared_ptr<ColumnDataWriter> mDebugWriter;
    boost::shared_ptr<ColumnDataWriter> mModeEventWriter;
    double mEventStartTime;
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
     * SetCellCycleDuration() method to set length of cell cycle (default 1.0 as Boije's model is purely generational)
     */
    void SetCellCycleDuration();

    /**
     * Overridden builder method to create new copies of
     * this cell-cycle model.
     *
     * @return new cell-cycle model
     */
    AbstractCellCycleModel* CreateCellCycleModel();

    /** Overridden ResetForDivision() method. */
    void ResetForDivision();

    /** Overridden InitialiseDaughterCell() method. */
    void InitialiseDaughterCell();

    void SetPostMitoticType(boost::shared_ptr<AbstractCellProperty> p_PostMitoticType);
    void SetModelParameters(const double normalMu, const double normalSigma, const double PP, const double PD);
    void SetModelProperties(boost::shared_ptr<AbstractCellProperty> p_RPh_Type,
                            boost::shared_ptr<AbstractCellProperty> p_AC_Type,
                            boost::shared_ptr<AbstractCellProperty> p_BC_Type,
                            boost::shared_ptr<AbstractCellProperty> p_MG_Type
                            );


    void EnableModeEventOutput(boost::shared_ptr<ColumnDataWriter> modeWriter, double eventStart, unsigned seed,
                               int timeID, std::vector<int> varIDs);
    void EnableModelDebugOutput(boost::shared_ptr<ColumnDataWriter> debugWriter);
    void EnableSequenceSampler(boost::shared_ptr<AbstractCellProperty> label);

    void WriteModeEventOutput();
    void WriteDebugData(double percentile);

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
