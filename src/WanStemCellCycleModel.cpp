#include "WanStemCellCycleModel.hpp"

WanStemCellCycleModel::WanStemCellCycleModel() :
        AbstractSimpleCellCycleModel(), mOutput(false), mEventStartTime(24.0), mDebug(false), mTimeID(), mVarIDs(), mDebugWriter(), mRPCfate(
                false), mGammaShift(4.0), mGammaShape(2.0), mGammaScale(1.0), mSeed(0), mTimeDependentCycleDuration(
                false), mPeakRateTime(), mIncreasingRateSlope(), mDecreasingRateSlope(), mBaseGammaScale(), mp_TransitType(), mHeParamVector(
                { 8, 15, 1, 0, .2, .4, .2, 0, 4, 2, 1, 1 })
{
}

WanStemCellCycleModel::WanStemCellCycleModel(const WanStemCellCycleModel& rModel) :
        AbstractSimpleCellCycleModel(rModel), mOutput(rModel.mOutput), mEventStartTime(rModel.mEventStartTime), mDebug(
                rModel.mDebug), mTimeID(rModel.mTimeID), mVarIDs(rModel.mVarIDs), mDebugWriter(rModel.mDebugWriter), mRPCfate(
                rModel.mRPCfate), mGammaShift(rModel.mGammaShift), mGammaShape(rModel.mGammaShape), mGammaScale(
                rModel.mGammaScale), mSeed(rModel.mSeed), mTimeDependentCycleDuration(
                rModel.mTimeDependentCycleDuration), mPeakRateTime(rModel.mPeakRateTime), mIncreasingRateSlope(
                rModel.mIncreasingRateSlope), mDecreasingRateSlope(rModel.mDecreasingRateSlope), mBaseGammaScale(
                rModel.mBaseGammaScale), mp_TransitType(rModel.mp_TransitType), mHeParamVector(rModel.mHeParamVector)
{
}

AbstractCellCycleModel* WanStemCellCycleModel::CreateCellCycleModel()
{
    return new WanStemCellCycleModel(*this);
}

void WanStemCellCycleModel::SetCellCycleDuration()
{
    RandomNumberGenerator* p_random_number_generator = RandomNumberGenerator::Instance();

    /**************************************
     * CELL CYCLE DURATION RANDOM VARIABLE
     *************************************/
    //Wan stem cell cycle length determined by the same formula as He RPCs
    if (!mTimeDependentCycleDuration) //Normal operation, cell cycle length stays constant
    {
        //He cell cycle length determined by shifted gamma distribution reflecting 4 hr refractory period followed by gamma pdf
        mCellCycleDuration = mGammaShift + p_random_number_generator->GammaRandomDeviate(mGammaShape, mGammaScale);
    }

    /****
     * Variable cycle length
     * Give -ve mIncreasingRateSlope and +ve mDecreasingRateSlope,
     * cell cycle length linearly declines (increasing rate), then increases, switching at mPeakRateTime
     ****/
    else
    {
        double currTime = SimulationTime::Instance()->GetTime();
        if (currTime <= mPeakRateTime)
        {
            mGammaScale = mBaseGammaScale * currTime * mIncreasingRateSlope;
        }
        if (currTime > mPeakRateTime)
        {
            mGammaScale = (mBaseGammaScale * mPeakRateTime * mIncreasingRateSlope)
                    + (mBaseGammaScale * (currTime - mPeakRateTime) * mDecreasingRateSlope);
        }
    }
}

void WanStemCellCycleModel::ResetForDivision()
{
    /********************************************
     * RPC-fated cells are given HeCellCycleModel
     ********************************************/

    if (mRPCfate)
    {
        double tiLOffset = -(SimulationTime::Instance()->GetTime());
        //Initialise a HeCellCycleModel and set it up with appropriate TiL values
        HeCellCycleModel* p_cycle_model = new HeCellCycleModel;
        p_cycle_model->SetModelParameters(tiLOffset, mHeParamVector[0], mHeParamVector[1], mHeParamVector[2],
                                          mHeParamVector[3], mHeParamVector[4], mHeParamVector[5], mHeParamVector[6],
                                          mHeParamVector[7], mHeParamVector[8], mHeParamVector[9], mHeParamVector[10],
                                          mHeParamVector[11]);
        p_cycle_model->Initialise();
        mpCell->SetCellCycleModel(p_cycle_model);
        mpCell->SetCellProliferativeType(mp_TransitType);
    }
    else
    {
        /****************
         * Write mitotic event to relevant files
         * *************/
        if (mDebug)
        {
            WriteDebugData();
        }

        if (mOutput)
        {
            WriteModeEventOutput();
        }

        AbstractSimpleCellCycleModel::ResetForDivision();
    }
}

void WanStemCellCycleModel::InitialiseDaughterCell()
{
    mRPCfate = 1;
}

void WanStemCellCycleModel::SetModelParameters(double gammaShift, double gammaShape, double gammaScale,
                                               std::vector<double> heParamVector)
{

    mGammaShift = gammaShift;
    mGammaShape = gammaShape;
    mGammaScale = gammaScale;
    mHeParamVector = heParamVector;
}

void WanStemCellCycleModel::SetTimeDependentCycleDuration(double peakRateTime, double increasingSlope,
                                                          double decreasingSlope)
{
    mTimeDependentCycleDuration = true;
    mPeakRateTime = peakRateTime;
    mIncreasingRateSlope = increasingSlope;
    mDecreasingRateSlope = decreasingSlope;
    mBaseGammaScale = mGammaScale;
}

void WanStemCellCycleModel::SetTransitType(boost::shared_ptr<AbstractCellProperty> p_TransitType)
{
    mp_TransitType = p_TransitType;
}

void WanStemCellCycleModel::EnableModeEventOutput(double eventStart, unsigned seed)
{
    mOutput = true;
    mEventStartTime = eventStart;
    mSeed = seed;
}

void WanStemCellCycleModel::WriteModeEventOutput()
{
    double currentTime = SimulationTime::Instance()->GetTime() + mEventStartTime;
    CellPtr currentCell = GetCell();
    double currentCellID = (double) currentCell->GetCellId();
    (*LogFile::Instance()) << currentTime << "\t" << mSeed << "\t" << currentCellID << "\t" << "4" << "\n"; // "4" indicates stem cell asymmetric division mode
}

void WanStemCellCycleModel::EnableModelDebugOutput(boost::shared_ptr<ColumnDataWriter> debugWriter)
{
    mDebug = true;
    mDebugWriter = debugWriter;

    mTimeID = mDebugWriter->DefineUnlimitedDimension("Time", "h");
    mVarIDs.push_back(mDebugWriter->DefineVariable("CellID", "No"));
    mVarIDs.push_back(mDebugWriter->DefineVariable("TiL", "h"));
    mVarIDs.push_back(mDebugWriter->DefineVariable("CycleDuration", "h"));
    mVarIDs.push_back(mDebugWriter->DefineVariable("Phase2Boundary", "h"));
    mVarIDs.push_back(mDebugWriter->DefineVariable("Phase3Boundary", "h"));
    mVarIDs.push_back(mDebugWriter->DefineVariable("Phase", "No"));
    mVarIDs.push_back(mDebugWriter->DefineVariable("Dieroll", "Percentile"));
    mVarIDs.push_back(mDebugWriter->DefineVariable("MitoticMode", "Mode"));

    mDebugWriter->EndDefineMode();
}

void WanStemCellCycleModel::WriteDebugData()
{
    double currentTime = SimulationTime::Instance()->GetTime();
    CellPtr currentCell = GetCell();
    double currentCellID = (double) currentCell->GetCellId();

    mDebugWriter->PutVariable(mTimeID, currentTime);
    mDebugWriter->PutVariable(mVarIDs[0], currentCellID);
    mDebugWriter->PutVariable(mVarIDs[1], 0);
    mDebugWriter->PutVariable(mVarIDs[2], mCellCycleDuration);
    mDebugWriter->PutVariable(mVarIDs[3], 0);
    mDebugWriter->PutVariable(mVarIDs[4], 0);
    mDebugWriter->PutVariable(mVarIDs[5], 0);
    mDebugWriter->PutVariable(mVarIDs[7], 4); //4 indicates stem cell asymmetric division
    mDebugWriter->AdvanceAlongUnlimitedDimension();
}

/******************
 * UNUSED FUNCTIONS (required for class nonvirtuality, do not remove)
 ******************/

double WanStemCellCycleModel::GetAverageTransitCellCycleTime()
{
    return (0.0);
}

double WanStemCellCycleModel::GetAverageStemCellCycleTime()
{
    return (0.0);
}

void WanStemCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<CellCycleDuration>" << mCellCycleDuration << "</CellCycleDuration>\n";

    // Call method on direct parent class
    AbstractSimpleCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(WanStemCellCycleModel)
