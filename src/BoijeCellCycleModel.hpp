#ifndef BOIJECELLCYCLEMODEL_HPP_
#define BOIJECELLCYCLEMODEL_HPP_

#include "AbstractSimpleCellCycleModel.hpp"
#include "BoijeRetinalNeuralFates.hpp"
#include "RandomNumberGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

/**
 * A stochastic cell-cycle model where cells divide with a stochastic cell cycle duration
 * with the length of the cell cycle drawn from a uniform distribution
 * on [mMinCellCycleDuration, mMaxCellCycleDuration].
 *
 * If the cell is differentiated, then the cell cycle duration is set to be infinite,
 * so that the cell will never divide.
 */
class BoijeCellCycleModel : public AbstractSimpleCellCycleModel
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

        SerializableSingleton<RandomNumberGenerator>* p_wrapper = RandomNumberGenerator::Instance()->GetSerializationWrapper();
        archive & p_wrapper;
        archive & mCellCycleDuration;
        archive & mGeneration;
    }

protected:

    /** The generation of this cell (cells with a StemCellProliferativeType have a generation of 0) */
    unsigned mGeneration;
    boost::shared_ptr<AbstractCellProperty> mp_PostMitoticType;
    boost::shared_ptr<AbstractCellProperty> mp_RGC_Type;
    boost::shared_ptr<AbstractCellProperty> mp_AC_HC_Type;
    boost::shared_ptr<AbstractCellProperty> mp_PR_BC_Type;
    double mprobAtoh7;
	double mprobPtf1a;
	double mprobng;
	bool mAtoh7Signal;
	bool mPtf1aSignal;
	bool mNgSignal;


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
    BoijeCellCycleModel(const BoijeCellCycleModel& rModel);

public:

    /**
     * Constructor - just a default, mBirthTime is set in the AbstractCellCycleModel class.
     */
    BoijeCellCycleModel();

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

    /**
     * Sets the cell's generation.
     *
     * @param generation the cell's generation
     **/
    void SetGeneration(unsigned generation);

    /**
     * @return the cell's generation.
     */
    unsigned GetGeneration() const;  
    
    void SetPostMitoticType(boost::shared_ptr<AbstractCellProperty> p_PostMitoticType);
    void SetRGCType(boost::shared_ptr<AbstractCellProperty> p_RGC_Type);    
    void SetACHCType(boost::shared_ptr<AbstractCellProperty> p_AC_HC_Type);    
    void SetPRBCType(boost::shared_ptr<AbstractCellProperty> p_PR_BC_Type);
    void SetModelParameters(double probAtoh7, double probPtf1a, double probng);
    
    /**
     * Overridden GetAverageTransitCellCycleTime() method.
     *
     * @return the average of mMinCellCycleDuration and mMaxCellCycleDuration
     */
    double GetAverageTransitCellCycleTime();

    /**
     * Overridden GetAverageStemCellCycleTime() method.
     *
     * @return the average of mMinCellCycleDuration and mMaxCellCycleDuration
     */
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
CHASTE_CLASS_EXPORT(BoijeCellCycleModel)

#endif /*BOIJECELLCYCLEMODEL_HPP_*/
