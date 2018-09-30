#include "BoijeCellCycleModel.hpp"

BoijeCellCycleModel::BoijeCellCycleModel() :
        AbstractSimpleCellCycleModel(), mOutput(false), mEventStartTime(), mSequenceSampler(false), mSeqSamplerLabelSister(
                false), mDebug(false), mTimeID(), mVarIDs(), mDebugWriter(), mGeneration(0), mPhase2gen(3), mPhase3gen(5), mprobAtoh7(0.32), mprobPtf1a(
                0.30), mprobng(0.80), mAtoh7Signal(false), mPtf1aSignal(false), mNgSignal(false), mMitoticMode(0), mSeed(
                0), mp_PostMitoticType(), mp_RGC_Type(), mp_AC_HC_Type(), mp_PR_BC_Type(), mp_label_Type()
{
}

BoijeCellCycleModel::BoijeCellCycleModel(const BoijeCellCycleModel& rModel) :
        AbstractSimpleCellCycleModel(rModel), mOutput(rModel.mOutput), mEventStartTime(rModel.mEventStartTime), mSequenceSampler(
                rModel.mSequenceSampler), mSeqSamplerLabelSister(rModel.mSeqSamplerLabelSister), mDebug(rModel.mDebug), mTimeID(
                rModel.mTimeID), mVarIDs(rModel.mVarIDs), mDebugWriter(rModel.mDebugWriter), mGeneration(
                rModel.mGeneration), mPhase2gen(rModel.mPhase2gen), mPhase3gen(rModel.mPhase3gen), mprobAtoh7(rModel.mprobAtoh7), mprobPtf1a(rModel.mprobPtf1a), mprobng(
                rModel.mprobng), mAtoh7Signal(rModel.mAtoh7Signal), mPtf1aSignal(rModel.mPtf1aSignal), mNgSignal(
                rModel.mNgSignal), mMitoticMode(rModel.mMitoticMode), mSeed(rModel.mSeed), mp_PostMitoticType(
                rModel.mp_PostMitoticType), mp_RGC_Type(rModel.mp_RGC_Type), mp_AC_HC_Type(rModel.mp_AC_HC_Type), mp_PR_BC_Type(
                rModel.mp_PR_BC_Type), mp_label_Type(rModel.mp_label_Type)
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
    mGeneration++; //increment generation counter, so that the first division produces generation "1"
    RandomNumberGenerator* p_random_number_generator = RandomNumberGenerator::Instance();

    mMitoticMode = 0; //0=PP;1=PD;2=DD

    /************************************************
     * TRANSCRIPTION FACTOR RANDOM VARIABLES & RULES
     * 1st phase: only PP divisions (symmetric proliferative) permitted
     * 2nd phase: occurs
     ************************************************/

    double atoh7RV, ptf1aRV, ngRV = 0;

    if (mGeneration > mPhase2gen && mGeneration < mPhase3gen) //if the cell is in the 2nd model phase, all signals have nonzero probabilities at each division
    {
        //RVs take values evenly distributed across 0-1
        atoh7RV = p_random_number_generator->ranf();
        ptf1aRV = p_random_number_generator->ranf();
        ngRV = p_random_number_generator->ranf();

        if (atoh7RV < mprobAtoh7)
        {
            mAtoh7Signal = true;
        }
        else
        {
            mAtoh7Signal = false;
        }
        if (ptf1aRV < mprobPtf1a)
        {
            mPtf1aSignal = true;
        }
        else
        {
            mPtf1aSignal = false;
        }
        if (ngRV < mprobng)
        {
            mNgSignal = true;
        }
        else
        {
            mNgSignal = false;
        }
    }

    /****************
     * Symmetric postmitotic specification rules
     * -(asymmetric postmitotic rules specified in InitialiseDaughterCell();)
     * *************/

    if (mGeneration > mPhase3gen) //if the cell is in the 3rd model phase, only ng signal has a nonzero probability
    {
        ngRV = p_random_number_generator->ranf();
        ; //roll a probability die for the ng signal
        if (mAtoh7Signal == true)
        {
            mAtoh7Signal = false;
            mMitoticMode = 1; //Atoh7 signal gives an asymmetrical PD divison
        }
        if (mPtf1aSignal == true)
        {
            mPtf1aSignal = false;
        }
        if (ngRV < mprobng)
        {
            mNgSignal = true;
        }
        else
        {
            mNgSignal = false;
        }
    }

    if (mPtf1aSignal == true && mAtoh7Signal == false) //Ptf1A alone gives a symmetrical postmitotic AC/HC division
    {
        mMitoticMode = 2;
        mpCell->SetCellProliferativeType(mp_PostMitoticType);
        mpCell->AddCellProperty(mp_AC_HC_Type);
    }

    if (mPtf1aSignal == false && mAtoh7Signal == false && mNgSignal == true) //ng alone gives a symmetrical postmitotic PR/BC division
    {
        mMitoticMode = 2;
        mpCell->SetCellProliferativeType(mp_PostMitoticType);
        mpCell->AddCellProperty(mp_PR_BC_Type);
    }

    /****************
     * Write mitotic event to file if appropriate
     * mDebug: many files: detailed per-lineage info switch; intended for TestHeInductionCountFixture
     * mOutput: 1 file: time, seed, cellID, mitotic mode, intended for TestHeMitoticModeRateFixture
     * *************/

    if (mDebug)
    {
        WriteDebugData(atoh7RV,ptf1aRV,ngRV);
    }

    if (mOutput)
    {
        WriteModeEventOutput();
    }

    AbstractSimpleCellCycleModel::ResetForDivision();

    /******************
     * SEQUENCE SAMPLER
     ******************/

    if (mpCell->HasCellProperty<CellLabel>())
    {
        (*LogFile::Instance()) << mMitoticMode;

        double labelDie = p_random_number_generator->ranf();
        if (labelDie <= .5)
        {
            mSeqSamplerLabelSister = true;
            mpCell->RemoveCellProperty<CellLabel>();
        }
        else
        {
            mSeqSamplerLabelSister = false;
        }
    }
}

void BoijeCellCycleModel::InitialiseDaughterCell()
{
    //Asymmetric specification rules

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

    /******************
     * SEQUENCE SAMPLER
     ******************/
    if (!mSeqSamplerLabelSister)
    {
        mpCell->RemoveCellProperty<CellLabel>();
    }
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
void BoijeCellCycleModel::SetSpecifiedTypes(boost::shared_ptr<AbstractCellProperty> p_RGC_Type, boost::shared_ptr<AbstractCellProperty> p_AC_HC_Type, boost::shared_ptr<AbstractCellProperty> p_PR_BC_Type)
{
    mp_RGC_Type = p_RGC_Type;
    mp_AC_HC_Type = p_AC_HC_Type;
    mp_PR_BC_Type = p_PR_BC_Type;
}

void BoijeCellCycleModel::SetModelParameters(unsigned phase2gen, unsigned phase3gen, double probAtoh7, double probPtf1a, double probng)
{
    mPhase2gen = phase2gen;
    mPhase3gen = phase3gen;
    mprobAtoh7 = probAtoh7;
    mprobPtf1a = probPtf1a;
    mprobng = probng;
}

void BoijeCellCycleModel::EnableModeEventOutput(double eventStart, unsigned seed)
{
    mOutput = true;
    mEventStartTime = eventStart;
    mSeed = seed;

}

void BoijeCellCycleModel::WriteModeEventOutput()
{
    double currentTime = SimulationTime::Instance()->GetTime() + mEventStartTime;
    CellPtr currentCell = GetCell();
    double currentCellID = (double) currentCell->GetCellId();
    (*LogFile::Instance()) << currentTime << "\t" << mSeed << "\t" << currentCellID << "\t" << mMitoticMode << "\n";
}

void BoijeCellCycleModel::EnableSequenceSampler(boost::shared_ptr<AbstractCellProperty> label)
{
    mSequenceSampler = true;
    mp_label_Type = label;
}

void BoijeCellCycleModel::EnableModelDebugOutput(boost::shared_ptr<ColumnDataWriter> debugWriter)
{
    mDebug = true;
    mDebugWriter = debugWriter;

    mTimeID = mDebugWriter->DefineUnlimitedDimension("Time", "h");
    mVarIDs.push_back(mDebugWriter->DefineVariable("CellID", "No."));
    mVarIDs.push_back(mDebugWriter->DefineVariable("Generation", "No."));
    mVarIDs.push_back(mDebugWriter->DefineVariable("MitoticMode", "0=PP;1=PD;2=DD"));
    mVarIDs.push_back(mDebugWriter->DefineVariable("atoh7Set", "Percentile"));
    mVarIDs.push_back(mDebugWriter->DefineVariable("atoh7RV", "Percentile"));
    mVarIDs.push_back(mDebugWriter->DefineVariable("ptf1aSet", "Percentile"));
    mVarIDs.push_back(mDebugWriter->DefineVariable("ptf1aRV", "Percentile"));
    mVarIDs.push_back(mDebugWriter->DefineVariable("ngSet", "Percentile"));
    mVarIDs.push_back(mDebugWriter->DefineVariable("ngRV", "Percentile"));

    mDebugWriter->EndDefineMode();
}

void BoijeCellCycleModel::WriteDebugData(double atoh7RV, double ptf1aRV, double ngRV)
{
    double currentTime = SimulationTime::Instance()->GetTime();
    CellPtr currentCell = GetCell();
    double currentCellID = (double) currentCell->GetCellId();

    mDebugWriter->PutVariable(mTimeID, currentTime);
    mDebugWriter->PutVariable(mVarIDs[0], currentCellID);
    mDebugWriter->PutVariable(mVarIDs[1], mGeneration);
    mDebugWriter->PutVariable(mVarIDs[2], mMitoticMode);
    mDebugWriter->PutVariable(mVarIDs[3], mprobAtoh7);
    mDebugWriter->PutVariable(mVarIDs[4], atoh7RV);
    mDebugWriter->PutVariable(mVarIDs[5], mprobPtf1a);
    mDebugWriter->PutVariable(mVarIDs[6], ptf1aRV);
    mDebugWriter->PutVariable(mVarIDs[7], mprobng);
    mDebugWriter->PutVariable(mVarIDs[8], ngRV);
    mDebugWriter->AdvanceAlongUnlimitedDimension();
}

/******************
 * UNUSED FUNCTIONS (required for class nonvirtuality, do not remove)
 ******************/

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
