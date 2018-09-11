#include "../../../../projects/ISP/src/BoijeCellCycleModel.hpp"

#include "DifferentiatedCellProliferativeType.hpp"
#include <RandomNumberGenerator.hpp>
#include "BoijeRetinalNeuralFates.hpp"
#include "NeuralTFs.hpp"

BoijeCellCycleModel::BoijeCellCycleModel()
    : AbstractSimpleCellCycleModel(),
    mGeneration(0),
    mp_PostMitoticType(),
    mp_RGC_Type(),
    mp_AC_HC_Type(),
    mp_PR_BC_Type(),
    mprobAtoh7(0.32),
    mprobPtf1a(0.30),
    mprobng(0.80),
    mAtoh7Signal(false),
    mPtf1aSignal(false),
    mNgSignal(false)
{
}

BoijeCellCycleModel::BoijeCellCycleModel(const BoijeCellCycleModel& rModel)
   : AbstractSimpleCellCycleModel(rModel),
   mGeneration(rModel.mGeneration),
   mp_PostMitoticType(rModel.mp_PostMitoticType),
   mp_RGC_Type(rModel.mp_RGC_Type),
   mp_AC_HC_Type(rModel.mp_AC_HC_Type),
   mp_PR_BC_Type(rModel.mp_PR_BC_Type),
   mprobAtoh7(rModel.mprobAtoh7),
   mprobPtf1a(rModel.mprobPtf1a),
   mprobng(rModel.mprobng),
   mAtoh7Signal(),
   mPtf1aSignal(),
   mNgSignal()
{
}

AbstractCellCycleModel* BoijeCellCycleModel::CreateCellCycleModel()
{
    return new BoijeCellCycleModel(*this);
}

void BoijeCellCycleModel::SetCellCycleDuration()
{

	if (mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
		{
			mCellCycleDuration = DBL_MAX;
		}
		else
		{
			mCellCycleDuration = 1.0; 
		}

}

void BoijeCellCycleModel::ResetForDivision()
{
    mGeneration++; //increment generation counter, so that the first division is "1" etc
    RandomNumberGenerator* p_random_number_generator = RandomNumberGenerator::Instance();
	p_random_number_generator->Reseed(rand());
    
    /****************
     * Transcription factor signal rules
     * *************/
    
    if (mGeneration > 3 && mGeneration < 5) //if the cell is in the 2nd model phase, all signals have nonzero probabilities at each division
    {
		
		double atoh7Die = p_random_number_generator->ranf(); //roll a probability die for each transcription factor system
		double ptf1aDie = p_random_number_generator->ranf();
		double ngDie = p_random_number_generator->ranf();
		
		if (atoh7Die < mprobAtoh7) {mAtoh7Signal = true;} else {mAtoh7Signal = false;}
		if (ptf1aDie < mprobPtf1a) {mPtf1aSignal = true;} else {mPtf1aSignal = false;}
		if (ngDie < mprobng) {mNgSignal = true;} else {mNgSignal = false;}
	}

    /****************
     * Symmetric postmitotic specification rules
     * -(asymmetric postmitotic rules specified in InitialiseDaughterCell();)
     * *************/

	if (mGeneration > 5) //if the cell is in the 3rd model phase, only ng signal has a nonzero probability
	{
		double ngDie = p_random_number_generator->ranf();; //roll a probability die for the ng signal
		if (mAtoh7Signal == true) {mAtoh7Signal = false;}
		if (mPtf1aSignal == true) {mPtf1aSignal = false;}
		if (ngDie < mprobng) {mNgSignal = true;} else {mNgSignal = false;}
	}
	
	if (mPtf1aSignal == true && mAtoh7Signal == false) //Ptf1A alone gives a symmetrical postmitotic AC/HC division
	{
		mpCell->SetCellProliferativeType(mp_PostMitoticType);
		mpCell->AddCellProperty(mp_AC_HC_Type);
	}
	
	if (mPtf1aSignal == false && mAtoh7Signal == false && mNgSignal == true) //ng alone gives a symmetrical postmitotic PR/BC division
	{
		mpCell->SetCellProliferativeType(mp_PostMitoticType);
		mpCell->AddCellProperty(mp_PR_BC_Type);	
	}
    
	p_random_number_generator->Destroy();

    AbstractSimpleCellCycleModel::ResetForDivision();
}

void BoijeCellCycleModel::InitialiseDaughterCell()
{
	
	if (mAtoh7Signal == true)
	{
		if (mPtf1aSignal == true)
		{
			mpCell->SetCellProliferativeType(mp_PostMitoticType);
			mpCell->AddCellProperty(mp_AC_HC_Type);
		}
		else
		{
			mpCell->SetCellProliferativeType(mp_PostMitoticType);
			mpCell->AddCellProperty(mp_RGC_Type);
		}
	}
	
	
	AbstractSimpleCellCycleModel::InitialiseDaughterCell();
}

void BoijeCellCycleModel::SetGeneration(unsigned generation)
{
    mGeneration = generation;
}

unsigned BoijeCellCycleModel::GetGeneration() const
{
    return mGeneration;
}

void BoijeCellCycleModel::SetPostMitoticType(boost::shared_ptr<AbstractCellProperty> p_PostMitoticType)
{
		mp_PostMitoticType = p_PostMitoticType;
}
void BoijeCellCycleModel::SetRGCType(boost::shared_ptr<AbstractCellProperty> p_RGC_Type)
{
		mp_RGC_Type = p_RGC_Type;
}

void BoijeCellCycleModel::SetACHCType(boost::shared_ptr<AbstractCellProperty> p_AC_HC_Type)
{
		mp_AC_HC_Type = p_AC_HC_Type;
}

void BoijeCellCycleModel::SetPRBCType(boost::shared_ptr<AbstractCellProperty> p_PR_BC_Type)
{
		mp_PR_BC_Type = p_PR_BC_Type;
}

void BoijeCellCycleModel::SetModelParameters(double probAtoh7, double probPtf1a, double probng)
{
	mprobAtoh7 = probAtoh7;
	mprobPtf1a = probPtf1a;
	mprobng = probng;
}

double BoijeCellCycleModel::GetAverageTransitCellCycleTime()
{
    return mCellCycleDuration;
}

double BoijeCellCycleModel::GetAverageStemCellCycleTime()
{
    return mCellCycleDuration;
}

void BoijeCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<CellCycleDuration>" << mCellCycleDuration << "</CellCycleDuration>\n";

    // Call method on direct parent class
    AbstractSimpleCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(BoijeCellCycleModel)
