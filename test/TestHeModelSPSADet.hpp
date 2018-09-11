#include <cxxtest/TestSuite.h>
#include "Exception.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "PetscSetupAndFinalize.hpp"

#include <string>
#include <armadillo>
#include <valarray>
#include <iostream>
#include <cmath>

#include "../../projects/ISP/src/HeCellCycleModel.hpp"
#include "../../projects/ISP/src/OffLatticeSimulationPropertyStop.hpp"
#include "AbstractCellBasedTestSuite.hpp"

#include "CellsGenerator.hpp"
#include "WildTypeCellMutationState.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "NodesOnlyMesh.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "SmartPointers.hpp"
#include "ColumnDataWriter.hpp"

using namespace arma;

class TestHeModelSPSADet : public AbstractCellBasedTestSuite
{
public:

    /*******************
     * SPSA COEFFICIENTS
     *******************/

    const double a = .0035; //chosen to produce reasonable initial step sizes
    const double c = .1; //small positive constant; using same seeds for each loss function eval
    const double A = 5; //~10% of expected # iterations
    double alpha = .602;
    double gamma = .101;

    /*************************
     * SPSA & AIC UTILITY VARS
     *************************/
    const unsigned maxIterations = 50;
    const unsigned numberPoints = 3000; //total number of comparison points between model output and empirical data (30 per timepoint)
    const unsigned numberPointsPerInduction = numberPoints / 3; //numberPoints must be evenly divisible by 3!
    const unsigned thetaSize = 4; //4 parameters to optimise
    vec delta = vec(thetaSize); //vector for bernoulli random variables to perturb theta

    /************************
     * SIMULATION PARAMETERS
     ************************/

    //Define start and end RNG seeds; determines:
    //#lineages per loss function run
    //unique sequence of RNG results for each lineage
    const unsigned startSeed = 0;
    unsigned endSeed = 999;
    unsigned totalSeeds = endSeed - startSeed + 1;

    /************************
     *GLOBAL MODEL PARAMETERS
     ************************/

    //Values defining different marker induction timepoints & relative start time of TiL counter
    const double earliestLineageStartTime = 23.0; //RPCs begin to enter "He model regime" nasally at 23hpf
    const double latestLineageStartTime = 39.0; //last temporal retinal RPC has entered "He model regime" at 39hpf
    const std::vector<double> inductionTimes = { 24.0, 32.0, 48.0 };
    const double endTime = 72.0;

    /**************************
     * SPECIFIC MODEL PARAMETERS - THETAHAT
     **************************/
    //DETERMINISTIC MITOTIC MODE
    //Phase boundary shift parameters
    const double phase1Offset = 6;
    const double phase2Shape = 3;
    const double phase2Scale = 4;
    const double phaseSisterShiftWidths = .8;
    const unsigned DModelParams = 11;
    vec thetaD = {  phase1Offset, phase2Shape, phase2Scale, phaseSisterShiftWidths, };
    vec shiftScaleVec = { 1, .5, .5, .2}; //scales ak gain sequences for shift & scale variables

    /****************************
     * HE ET AL EMPIRICAL RESULTS
     ****************************/

    const vec he24counts = { 0, 0, 1, 0, 1, 2, 0, 7, 4, 9, 6, 5, 5, 6, 5, 5, 1, 2, 2, 2, 0, 1 };
    const vec he32counts = { 6, 20, 25, 22, 24, 17, 11, 7, 8, 5, 6, 3, 1, 6, 4, 1, 0, 0, 0, 2, 0, 0, 0, 1 };
    const vec he48counts = { 59, 86, 2, 12, 1, 2, 0, 1 };

    const unsigned lineagesSampled24 = 64;
    const unsigned lineagesSampled32 = 169;
    const unsigned lineagesSampled48 = 163;

    vec he24AICvec = vec(numberPointsPerInduction, fill::zeros);
    vec he32AICvec = vec(numberPointsPerInduction, fill::zeros);
    vec he48AICvec = vec(numberPointsPerInduction, fill::zeros);
    uvec centers = uvec(numberPointsPerInduction, fill::zeros);

    /*********************
     * TEST FIXTURES
     *********************/

    void TestHeModelSPSADetFixture()
    {
        Timer::Reset();
        Timer::Print("Begin TestHeSPSADet Fixture @ ");

        //Setup params logfile
        LogFile* p_log = LogFile::Instance();
        p_log->Set(0, "HeModelParameterOptimise", "DetParams");
        *p_log << "k\tPhase1Offset\tPhase2Shape\tPhase2Scale\tphaseSisterShiftWidths\n";

        //Set up AIC comparison vectors
        for (unsigned i = 0; i < he24counts.n_elem; i++)
        {
            he24AICvec(i) = he24counts(i);
        }
        for (unsigned i = 0; i < he32counts.n_elem; i++)
        {
            he32AICvec(i) = he32counts(i);
        }
        for (unsigned i = 0; i < he48counts.n_elem; i++)
        {
            he48AICvec(i) = he48counts(i);
        }
        for (unsigned i = 0; i < centers.n_elem; i++)
        {
            centers(i) = i + 1;
        }

        unsigned k = 0;

        while (k < maxIterations)
        {

            //at 35th iterate, increase seeds to 5000
            if (k == 35)
            {
                endSeed = 4999;
                totalSeeds = endSeed + 1;
            }

            *p_log << k << "\t" << thetaD(0) << "\t" << thetaD(1) << "\t" << thetaD(2) << "\t"
                    << thetaD(3) << "\n";

            getBernoulliVector(); //populate deltak perturbation vector

            double ak = a / (pow((A + k + 1), alpha)); //calculate ak from gain sequence
            vec scaledak = ak * shiftScaleVec; // scale ak appropriately for gain & shift params

            double ck = c / (pow((k + 1), gamma)); //calculate ck from gain sequence

            vec projectedTheta = project(thetaD, ck);

            //Calculate theta+ and theta- vectors for gradient estimate
            vec thetaPlus = projectedTheta + ck * delta;
            vec thetaMinus = projectedTheta - ck * delta;

            double positiveAIC = evaluateAIC(thetaPlus, true, DModelParams, numberPoints);
            double negativeAIC = evaluateAIC(thetaMinus, true, DModelParams, numberPoints);

            vec ghat = ((positiveAIC - negativeAIC) / (2 * ck)) * delta;

            *p_log << "PositiveAIC: " << positiveAIC << " NegativeAIC: " << negativeAIC << "\n";
            *p_log << "ghat0: " << ghat(0) << "\n";

            thetaD = thetaD - scaledak % ghat;

            //project adjusted thetaD
            thetaD = project(thetaD, 0);

            k++;
        }

        *p_log << "Final thetaD:\n" << k << "\t" << thetaD(0) << "\t" << thetaD(1) << "\t" << thetaD(2) << "\t"
                << thetaD(3) << "\n";

        p_log->Close();
    }

    /*******************
     * UTILITY FUNCTIONS
     *******************/

    void getBernoulliVector()
    {
        RandomNumberGenerator* p_RNG = RandomNumberGenerator::Instance();
        p_RNG->Reseed(rand());

        for (unsigned i = 0; i < thetaSize; i++)
        {
            double randVar = p_RNG->ranf();
            if (randVar < .5) delta(i) = -1;
            else delta(i) = 1;
        }
    }

    /******************************
     * SPSA THETA CONSTRAINTS
     ******************************/

    vec project(vec theta, double boundary)
    {
        for (unsigned i = 0; i < thetaSize; i++)
        {
            if (theta[i] < 0.0001 + boundary) theta[i] = 0.0001 + boundary;
        }

        return theta;
    }

    double evaluateAIC(vec theta, bool detMode, unsigned numberParams, unsigned totalPoints)
    {

        //vector to store the 3 residual sums of squares
        vec RSS = vec(3);
        unsigned RSSindex = 0;
        totalSeeds = endSeed + 1;

        //iterate through inductionTimes
        for (auto inductionTime : inductionTimes)
        {
            vec empiricalCounts;
            unsigned empiricalLineagesSampled;

            if (inductionTime == 24)
            {
                empiricalCounts = he24AICvec;
                empiricalLineagesSampled = lineagesSampled24;
            }
            if (inductionTime == 32)
            {
                empiricalCounts = he32AICvec;
                empiricalLineagesSampled = lineagesSampled32;
            }
            if (inductionTime == 48)
            {
                empiricalCounts = he48AICvec;
                empiricalLineagesSampled = lineagesSampled48;
            }

            uvec histo = hist(runSimulator(inductionTime, theta, detMode), centers);
            vec probHisto = conv_to<vec>::from(histo);
            probHisto = probHisto / totalSeeds;
            vec probEmpirical = empiricalCounts / empiricalLineagesSampled;
            vec residual = probHisto - probEmpirical;
            vec residualSquares = square(residual);
            RSS(RSSindex) = sum(residualSquares);

            RSSindex++;
        }

        double AIC = 2 * numberParams + totalPoints * log(sum(RSS));
        return AIC;
    }

    /**************************************
     * SIMULATOR- returns vector of counts
     **************************************/

    uvec runSimulator(double inductionTime, vec theta, bool deterministic)
    {
        uvec counts = uvec(totalSeeds);

        RandomNumberGenerator* p_RNG = RandomNumberGenerator::Instance();

        //iterate through supplied seed range, executing one simulation per seed
        for (unsigned seed = startSeed; seed <= endSeed; seed++)
        {

            //Initialise pointers to relevant ProliferativeTypes and Properties
            MAKE_PTR(WildTypeCellMutationState, p_state);
            MAKE_PTR(TransitCellProliferativeType, p_Mitotic);
            MAKE_PTR(DifferentiatedCellProliferativeType, p_PostMitotic);

            //Reseed the RNG with the required seed
            p_RNG->Reseed(seed);

            //Initialise a HeCellCycleModel and set it up with appropriate TiL values, also obtaining the correct simulation end time
            HeCellCycleModel* p_cycle_model = new HeCellCycleModel;
            double simEndTime = SetupNTHeModel(p_cycle_model, inductionTime, theta, p_PostMitotic, deterministic);

            //Setup vector containing lineage founder*/
            std::vector<CellPtr> cells;
            CellPtr p_cell(new Cell(p_state, p_cycle_model));
            p_cell->SetCellProliferativeType(p_Mitotic);
            p_cell->InitialiseCellCycleModel();
            cells.push_back(p_cell);

            //Generate 1x1 mesh for single-cell colony
            HoneycombMeshGenerator generator(1, 1);
            MutableMesh<2, 2>* p_generating_mesh = generator.GetMesh();
            //Prepare nodes-only mesh
            NodesOnlyMesh<2> mesh;
            mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);

            //Setup cell population
            NodeBasedCellPopulation<2>* cell_population(new NodeBasedCellPopulation<2>(mesh, cells));

            //Setup simulator & run simulation
            boost::shared_ptr<OffLatticeSimulationPropertyStop<2>> p_simulator(
                    new OffLatticeSimulationPropertyStop<2>(*cell_population));
            p_simulator->SetStopProperty(p_Mitotic); //simulation to stop if no mitotic cells are left
            p_simulator->SetDt(0.05);
            p_simulator->SetEndTime(simEndTime);
            p_simulator->SetOutputDirectory("UnusedSimOutputOptimisationsDet"); //unused output
            p_simulator->Solve();

            //Count lineage size
            unsigned count = cell_population->GetNumRealCells();
            counts(seed) = count;
            double currSimTime = SimulationTime::Instance()->GetTime();

            //Output simulation report to console
            std::stringstream reportStream;
            reportStream << "inductiontime " << inductionTime << "Simulation Seed " << seed << ": Final Lineage Count: "
                    << count << " Sim Complete; Last SimTime " << currSimTime << " @ ";
            auto simReport = reportStream.str();
            Timer::Print(simReport);

            //Reset for next simulation
            SimulationTime::Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);
            delete cell_population;
        }

        p_RNG->Destroy();

        return counts;
    }

    double SetupNTHeModel(HeCellCycleModel* p_model, double inductionTime, vec theta,
                          boost::shared_ptr<AbstractCellProperty> pmType, bool deterministic = false,
                          std::string debugDir = "debug")
    {
        //get the RNG instance
        RandomNumberGenerator* p_RNG = RandomNumberGenerator::Instance();
        /******************************************************
         * Generate Time in Lineage figure from even distribution across nasal-temporal gradient
         ******************************************************/
        //generate random lineage start time from even random distro across earliest-latest start time figures
        double lineageStartTime = (p_RNG->ranf() * (latestLineageStartTime - earliestLineageStartTime))
                + earliestLineageStartTime;
        double currTiL;
        double currSimEndTime;
        //if the lineage started before the present induction time, give the cell positive TiL and run the full simulation time
        //this reflects induction of cells after the lineages' first mitosis
        if (lineageStartTime < inductionTime)
        {
            currTiL = inductionTime - lineageStartTime;
            currSimEndTime = endTime - inductionTime;
        }
        //if the lineage starts after the induction time, give it zero TiL run the appropriate-length simulation
        //(ie. the endTime is reduced by the amount of time after induction that the first mitosis occurs)
        if (lineageStartTime >= inductionTime)
        {
            currTiL = 0.0;
            currSimEndTime = endTime - lineageStartTime;
        }

        //Setup lineages' cycle model with appropriate parameters
        p_model->SetDimension(2);
        p_model->SetPostMitoticType(pmType);

        if (!deterministic)
        {
            /*p_model->SetModelParameters(currTiL, theta(0), theta(0) + theta(1),
                                        theta(2), phase1PD, theta(3), theta(4), theta(5), phase3PD);*/
        }
        else
        {
            //Gamma-distribute phase3 boundary
            double currPhase2Boundary = theta(0);
            double currPhase3Boundary = currPhase2Boundary + p_RNG->GammaRandomDeviate(theta(1), theta(2));

            p_model->SetDeterministicMode(currTiL, currPhase2Boundary,
                                          currPhase3Boundary, theta(3));
        }
        return currSimEndTime;
    }

}
;
