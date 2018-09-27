#include "GomesCellCycleModel.hpp"

GomesCellCycleModel::GomesCellCycleModel() :
        AbstractSimpleCellCycleModel(), mDebug(false), mOutput(false), mSequenceSampler(false), mSeqSamplerLabelSister(
                false), mTimeID(), mVarIDs(), mDebugWriter(), mModeEventWriter(), mEventStartTime(), mNormalMu(3.9716), mNormalSigma(
                0.32839), mPP(.055), mPD(0.221), mpBC(.128), mpAC(.106), mpMG(.028), mMitoticMode(), mSeed(), mp_PostMitoticType(), mp_RPh_Type(), mp_BC_Type(), mp_AC_Type(), mp_MG_Type(), mp_label_Type()
{
}

GomesCellCycleModel::GomesCellCycleModel(const GomesCellCycleModel& rModel) :
        AbstractSimpleCellCycleModel(rModel), mDebug(rModel.mDebug), mOutput(rModel.mOutput), mSequenceSampler(
                rModel.mSequenceSampler), mSeqSamplerLabelSister(rModel.mSeqSamplerLabelSister), mTimeID(
                rModel.mTimeID), mVarIDs(rModel.mVarIDs), mDebugWriter(rModel.mDebugWriter), mModeEventWriter(
                rModel.mModeEventWriter), mEventStartTime(rModel.mEventStartTime), mNormalMu(rModel.mNormalMu), mNormalSigma(
                rModel.mNormalSigma), mPP(rModel.mPP), mPD(rModel.mPD), mpBC(rModel.mpBC), mpAC(rModel.mpAC), mpMG(
                rModel.mpMG), mMitoticMode(rModel.mMitoticMode), mSeed(rModel.mSeed), mp_PostMitoticType(
                rModel.mp_PostMitoticType), mp_RPh_Type(rModel.mp_RPh_Type), mp_BC_Type(rModel.mp_BC_Type), mp_AC_Type(
                rModel.mp_AC_Type), mp_MG_Type(rModel.mp_MG_Type), mp_label_Type(rModel.mp_label_Type)
{
}

AbstractCellCycleModel* GomesCellCycleModel::CreateCellCycleModel()
{
    return new GomesCellCycleModel(*this);
}

void GomesCellCycleModel::SetCellCycleDuration()
{
    RandomNumberGenerator* p_random_number_generator = RandomNumberGenerator::Instance();

    //Gomes cell cycle length determined by lognormal distribution with default mean 56 hr, std 18.9 hrs.
    mCellCycleDuration = exp(p_random_number_generator->NormalRandomDeviate(mNormalMu, mNormalSigma));
}

void GomesCellCycleModel::ResetForDivision()
{
    /****************
     * Mitotic mode rules
     * *************/
    RandomNumberGenerator* p_random_number_generator = RandomNumberGenerator::Instance();

    //Check time in lineage and determine current mitotic mode phase
    mMitoticMode = 0; //0=PP;1=PD;2=DD

    /******************************
     * MITOTIC MODE RANDOM VARIABLE
     ******************************/
    //initialise mitoticmode random variable, set mitotic mode appropriately after comparing to mode probability array
    double mitoticModeDie = p_random_number_generator->ranf();

    if (mitoticModeDie > mPP && mitoticModeDie <= mPP + mPD)
    {
        mMitoticMode = 1;
    }
    if (mitoticModeDie > mPP + mPD)
    {
        mMitoticMode = 2;
    }

    /****************
     * Write mitotic event to file if appropriate
     * mDebug: many files: detailed per-lineage info switch; intended for TestHeInductionCountFixture
     * mOutput: 1 file: time, seed, cellID, mitotic mode, intended for TestHeMitoticModeRateFixture
     * *************/

    if (mDebug)
    {
        WriteDebugData(mitoticModeDie);
    }

    if (mOutput)
    {
        WriteModeEventOutput();
    }

    //set new cell cycle length (will be overwritten with DBL_MAX for DD divisions)
    AbstractSimpleCellCycleModel::ResetForDivision();

    /****************
     * Symmetric postmitotic specification rule
     * -(asymmetric postmitotic rule specified in InitialiseDaughterCell();)
     * *************/

    if (mMitoticMode == 2)
    {
        mpCell->SetCellProliferativeType(mp_PostMitoticType);
        mCellCycleDuration = DBL_MAX;
        /*********************
         * SPECIFICATION RULES
         ********************/
        double specificationDie = p_random_number_generator->ranf();
        if (specificationDie <= mpMG)
        {
            mpCell->AddCellProperty(mp_MG_Type);
        }
        if (specificationDie > mpMG && specificationDie <= mpMG + mpAC)
        {
            mpCell->AddCellProperty(mp_AC_Type);
        }
        if (specificationDie > mpMG + mpAC && specificationDie <= mpMG + mpAC + mpBC)
        {
            mpCell->AddCellProperty(mp_BC_Type);
        }
        if (specificationDie > mpMG + mpAC + mpBC)
        {
            mpCell->AddCellProperty(mp_RPh_Type);
        }
    }

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

void GomesCellCycleModel::InitialiseDaughterCell()
{

    //daughter cell's mCellCycleDuration is copied from parent; reset to new value from gamma PDF here
    SetCellCycleDuration();

    /************
     * PD-type division
     **********/


    if (mMitoticMode == 1)
    {
        RandomNumberGenerator* p_random_number_generator = RandomNumberGenerator::Instance();
        mpCell->SetCellProliferativeType(mp_PostMitoticType);
        mCellCycleDuration = DBL_MAX;
        /*********************
         * SPECIFICATION RULES
         ********************/
        double specificationDie = p_random_number_generator->ranf();
        if (specificationDie <= mpMG)
        {
            mpCell->AddCellProperty(mp_MG_Type);
        }
        if (specificationDie > mpMG && specificationDie <= mpMG + mpAC)
        {
            mpCell->AddCellProperty(mp_AC_Type);
        }
        if (specificationDie > mpMG + mpAC && specificationDie <= mpMG + mpAC + mpBC)
        {
            mpCell->AddCellProperty(mp_BC_Type);
        }
        if (specificationDie > mpMG + mpAC + mpBC)
        {
            mpCell->AddCellProperty(mp_RPh_Type);
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

void GomesCellCycleModel::SetModelParameters(const double normalMu, const double normalSigma, const double PP,
                                             const double PD, const double pBC, const double pAC, const double pMG)
{
    mNormalMu = normalMu;
    mNormalSigma = normalSigma;
    mPP = PP;
    mPD = PD;
    mpBC = pBC;
    mpAC = pAC;
    mpMG = pMG;

}

void GomesCellCycleModel::SetModelProperties(boost::shared_ptr<AbstractCellProperty> p_RPh_Type,
                                             boost::shared_ptr<AbstractCellProperty> p_AC_Type,
                                             boost::shared_ptr<AbstractCellProperty> p_BC_Type,
                                             boost::shared_ptr<AbstractCellProperty> p_MG_Type)
{
    mp_RPh_Type = p_RPh_Type;
    mp_AC_Type = p_AC_Type;
    mp_BC_Type = p_BC_Type;
    mp_MG_Type = p_MG_Type;
}

void GomesCellCycleModel::SetPostMitoticType(boost::shared_ptr<AbstractCellProperty> p_PostMitoticType)
{
    mp_PostMitoticType = p_PostMitoticType;
}

void GomesCellCycleModel::EnableModeEventOutput(double eventStart, unsigned seed)
{
    mOutput = true;
    mEventStartTime = eventStart;
    mSeed = seed;

}

void GomesCellCycleModel::WriteModeEventOutput()
{
    double currentTime = SimulationTime::Instance()->GetTime() + mEventStartTime;
    CellPtr currentCell = GetCell();
    double currentCellID = (double) currentCell->GetCellId();
    (*LogFile::Instance()) << currentTime << "\t" << mSeed << "\t" << currentCellID << "\t" << mMitoticMode << "\n";
}

void GomesCellCycleModel::EnableSequenceSampler(boost::shared_ptr<AbstractCellProperty> label)
{
    mSequenceSampler = true;
    mp_label_Type = label;
}


void GomesCellCycleModel::EnableModelDebugOutput(boost::shared_ptr<ColumnDataWriter> debugWriter)
{
    mDebug = true;
    mDebugWriter = debugWriter;

    mTimeID = mDebugWriter->DefineUnlimitedDimension("Time", "h");
    mVarIDs.push_back(mDebugWriter->DefineVariable("CellID", "No"));
    mVarIDs.push_back(mDebugWriter->DefineVariable("CycleDuration", "h"));
    mVarIDs.push_back(mDebugWriter->DefineVariable("PP", "Percentile"));
    mVarIDs.push_back(mDebugWriter->DefineVariable("PD", "Percentile"));
    mVarIDs.push_back(mDebugWriter->DefineVariable("Dieroll", "Percentile"));
    mVarIDs.push_back(mDebugWriter->DefineVariable("MitoticMode", "Mode"));

    mDebugWriter->EndDefineMode();
}

void GomesCellCycleModel::WriteDebugData(double percentileRoll)
{
    double currentTime = SimulationTime::Instance()->GetTime();
    CellPtr currentCell = GetCell();
    double currentCellID = (double) currentCell->GetCellId();

    mDebugWriter->PutVariable(mTimeID, currentTime);
    mDebugWriter->PutVariable(mVarIDs[0], currentCellID);
    mDebugWriter->PutVariable(mVarIDs[1], mCellCycleDuration);
    mDebugWriter->PutVariable(mVarIDs[2], mPP);
    mDebugWriter->PutVariable(mVarIDs[3], mPD);
    mDebugWriter->PutVariable(mVarIDs[4], percentileRoll);
    mDebugWriter->PutVariable(mVarIDs[5], mMitoticMode);
    mDebugWriter->AdvanceAlongUnlimitedDimension();
}

/******************
 * UNUSED FUNCTIONS (required for class nonvirtuality, do not remove)
 ******************/

double GomesCellCycleModel::GetAverageTransitCellCycleTime()
{
    return (0.0);
}

double GomesCellCycleModel::GetAverageStemCellCycleTime()
{
    return (0.0);
}

void GomesCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<CellCycleDuration>" << mCellCycleDuration << "</CellCycleDuration>\n";

    // Call method on direct parent class
    AbstractSimpleCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(GomesCellCycleModel)
