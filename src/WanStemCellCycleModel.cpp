#include "WanStemCellCycleModel.hpp"

WanStemCellCycleModel::WanStemCellCycleModel() :
        AbstractSimpleCellCycleModel(), mExpandingStemPopulation(false), mPopulation(), mOutput(false), mEventStartTime(
                72.0), mDebug(false), mTimeID(), mVarIDs(), mDebugWriter(), mBasePopulation(), mGammaShift(4.0), mGammaShape(
                2.0), mGammaScale(1.0), mMitoticMode(0), mSeed(0), mTimeDependentCycleDuration(false), mPeakRateTime(), mIncreasingRateSlope(), mDecreasingRateSlope(), mBaseGammaScale(), mHeParamVector(
                { 8, 15, 1, 0, .2, .4, .2, 0, 4, 2, 1, 1 })
{
}

WanStemCellCycleModel::WanStemCellCycleModel(const WanStemCellCycleModel& rModel) :
        AbstractSimpleCellCycleModel(rModel), mExpandingStemPopulation(rModel.mExpandingStemPopulation), mPopulation(
                rModel.mPopulation), mOutput(rModel.mOutput), mEventStartTime(rModel.mEventStartTime), mDebug(
                rModel.mDebug), mTimeID(rModel.mTimeID), mVarIDs(rModel.mVarIDs), mDebugWriter(rModel.mDebugWriter), mBasePopulation(
                rModel.mBasePopulation), mGammaShift(rModel.mGammaShift), mGammaShape(rModel.mGammaShape), mGammaScale(
                rModel.mGammaScale), mMitoticMode(rModel.mMitoticMode), mSeed(rModel.mSeed), mTimeDependentCycleDuration(
                rModel.mTimeDependentCycleDuration), mPeakRateTime(rModel.mPeakRateTime), mIncreasingRateSlope(
                rModel.mIncreasingRateSlope), mDecreasingRateSlope(rModel.mDecreasingRateSlope), mBaseGammaScale(
                rModel.mBaseGammaScale), mHeParamVector(rModel.mHeParamVector)
{
}

AbstractCellCycleModel* WanStemCellCycleModel::CreateCellCycleModel()
{
    return new WanStemCellCycleModel(*this);
}

void WanStemCellCycleModel::SetCellCycleDuration()
{
    RandomNumberGenerator* p_random_number_generator = RandomNumberGenerator::Instance();

    mCellCycleDuration = mGammaShift + p_random_number_generator->GammaRandomDeviate(mGammaShape, mGammaScale);

    /**************************************
     * CELL CYCLE DURATION RANDOM VARIABLE
     *************************************/
    //Wan stem cell cycle length determined by the same formula as He RPCs
    /*
     if (!mTimeDependentCycleDuration) //Normal operation, cell cycle length stays constant
     {
     //He cell cycle length determined by shifted gamma distribution reflecting 4 hr refractory period followed by gamma pdf
     mCellCycleDuration = mGammaShift + p_random_number_generator->GammaRandomDeviate(mGammaShape, mGammaScale);
     }*/

    /****
     * Variable cycle length
     * Give -ve mIncreasingRateSlope and +ve mDecreasingRateSlope,
     * cell cycle length linearly declines (increasing rate), then increases, switching at mPeakRateTime
     ****/
    /*
     else
     {
     double currTime = SimulationTime::Instance()->GetTime();
     if (currTime <= mPeakRateTime)
     {
     mGammaScale = std::max((mBaseGammaScale - currTime * mIncreasingRateSlope), .0000000000001);
     }
     if (currTime > mPeakRateTime)
     {
     mGammaScale = std::max(((mBaseGammaScale - mPeakRateTime * mIncreasingRateSlope)
     + (mBaseGammaScale + (currTime - mPeakRateTime) * mDecreasingRateSlope)),.0000000000001);
     }
     Timer::Print("mGammaScale: " + std::to_string(mGammaScale));
     mCellCycleDuration = mGammaShift + p_random_number_generator->GammaRandomDeviate(mGammaShape, mGammaScale);
     }
     */
}

void WanStemCellCycleModel::ResetForDivision()
{

    mMitoticMode = 1; //by default, asymmetric division giving rise to He cell (mode 1)

    if (mExpandingStemPopulation)
    {
        double currRetinaAge = SimulationTime::Instance()->GetTime() + mEventStartTime;
        double lensGrowthFactor = .09256 * pow(currRetinaAge, .52728); // power law model fit for lens growth
        unsigned currentPopulationTarget = int(std::round(mBasePopulation * lensGrowthFactor));

        unsigned currentStemPopulation = (mPopulation->GetCellProliferativeTypeCount())[0];
        if (currentStemPopulation < currentPopulationTarget)
        {
            mMitoticMode = 0; //if the current population is < target, symmetrical stem-stem division occurs (mode 0)
        }
    }

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

void WanStemCellCycleModel::Initialise()
{

    boost::shared_ptr<AbstractCellProperty> p_Stem =
            mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<StemCellProliferativeType>();
    mpCell->SetCellProliferativeType(p_Stem);

    SetCellCycleDuration();
}

void WanStemCellCycleModel::InitialiseDaughterCell()
{

    if (mMitoticMode == 1)
    {
        /********************************************
         * RPC-fated cells are given HeCellCycleModel
         ********************************************/

        double tiLOffset = -(SimulationTime::Instance()->GetTime());
        //Initialise a HeCellCycleModel and set it up with appropriate TiL value & parameters
        HeCellCycleModel* p_cycle_model = new HeCellCycleModel;
        p_cycle_model->SetModelParameters(tiLOffset, mHeParamVector[0], mHeParamVector[1], mHeParamVector[2],
                                          mHeParamVector[3], mHeParamVector[4], mHeParamVector[5], mHeParamVector[6],
                                          mHeParamVector[7], mHeParamVector[8], mHeParamVector[9], mHeParamVector[10],
                                          mHeParamVector[11]);
        p_cycle_model->EnableKillSpecified();

        //if debug output is enabled for the stem cell, enable it for its progenitor offspring
        if (mDebug)
        {
            p_cycle_model->PassDebugWriter(mDebugWriter, mTimeID, mVarIDs);
        }

        mpCell->SetCellCycleModel(p_cycle_model);
        p_cycle_model->Initialise();
    }

    else
    {
        SetCellCycleDuration();
    }
}

void WanStemCellCycleModel::SetModelParameters(double gammaShift, double gammaShape, double gammaScale,
                                               std::vector<double> heParamVector)
{
    mGammaShift = gammaShift;
    mGammaShape = gammaShape;
    mGammaScale = gammaScale;
    mHeParamVector = heParamVector;
}

void WanStemCellCycleModel::EnableExpandingStemPopulation(int basePopulation,
                                                          boost::shared_ptr<AbstractCellPopulation<2>> p_population)
{
    mExpandingStemPopulation = true;
    mBasePopulation = basePopulation;
    mPopulation = p_population;
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
    (*LogFile::Instance()) << currentTime << "\t" << mSeed << "\t" << currentCellID << "\t" << mMitoticMode << "\n"; //
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
    mVarIDs.push_back(mDebugWriter->DefineVariable("MitoticModeRV", "Percentile"));
    mVarIDs.push_back(mDebugWriter->DefineVariable("MitoticMode", "Mode"));
    mVarIDs.push_back(mDebugWriter->DefineVariable("Label", "binary"));

    mDebugWriter->EndDefineMode();
}

void WanStemCellCycleModel::WriteDebugData()
{
    double currentTime = SimulationTime::Instance()->GetTime();
    double currentCellID = mpCell->GetCellId();

    mDebugWriter->PutVariable(mTimeID, currentTime);
    mDebugWriter->PutVariable(mVarIDs[0], currentCellID);
    mDebugWriter->PutVariable(mVarIDs[1], 0);
    mDebugWriter->PutVariable(mVarIDs[2], mCellCycleDuration);
    mDebugWriter->PutVariable(mVarIDs[3], 0);
    mDebugWriter->PutVariable(mVarIDs[4], 0);
    mDebugWriter->PutVariable(mVarIDs[5], 0);
    mDebugWriter->PutVariable(mVarIDs[7], mMitoticMode);
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
