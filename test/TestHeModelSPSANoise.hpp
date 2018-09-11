#include <cxxtest/TestSuite.h>
#include "Exception.hpp"
#include <string>
#include <armadillo>
#include <valarray>
#include <iostream>
#include <cmath>

#include "../../projects/ISP/src/HeCellCycleModel.hpp"
#include "../../projects/ISP/src/OffLatticeSimulationPropertyStop.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "PetscSetupAndFinalize.hpp"

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

class TestHeModelSPSA : public AbstractCellBasedTestSuite
{
public:

    /*************************
     * AIC UTILITY VARS
     *************************/
    const unsigned maxIterations = 5;
    const unsigned numberPoints = 3000; //total number of comparison points between model output and empirical data (30 per timepoint)
    const unsigned numberPointsPerInduction = numberPoints / 3; //numberPoints must be evenly divisible by 3!
    vec AICvec = vec(maxIterations);

    /************************
     * SIMULATION PARAMETERS
     ************************/

    unsigned totalSeeds = 1000;

    /************************
     *GLOBAL MODEL PARAMETERS
     ************************/

    //Values defining different marker induction timepoints & relative start time of TiL counter
    const double earliestLineageStartTime = 23.0; //RPCs begin to enter "He model regime" nasally at 23hpf
    const double latestLineageStartTime = 39.0; //last temporal retinal RPC has entered "He model regime" at 39hpf
    const std::vector<double> inductionTimes = { 24.0, 32.0, 48.0 };
    const double endTime = 72.0;
    //Cell cycle length gamma PDF parameters
    const double gammaShift = 4.0;
    const double gammaShape = 2.0;
    const double gammaScale = 1.0;

    /**************************
     * SPECIFIC MODEL PARAMETERS - THETAHAT
     **************************/
    //STOCHASTIC MITOTIC MODE
    //mitotic mode per-phase probabilities
    //3-phase mitotic mode time periodisation
    const double mitoticModePhase2 = 8.0;
    const double mitoticModePhase3 = 7.0;
    const double phase1PP = 1.0;
    const double phase1PD = 0.0;
    const double phase2PP = 0.2;
    const double phase2PD = 0.4;
    const double phase3PP = 0.2;
    const double phase3PD = 0.0;
    const unsigned HeModelParams = 15;
    vec thetaS = { mitoticModePhase2, mitoticModePhase3, phase1PP, phase2PP, phase2PD, phase3PP };
    //DETERMINISTIC MITOTIC MODE
    //Phase boundary shift parameters
    const double phaseSisterShiftWidths = 0.8;
    const double b1Shape = 2.0;
    const double b1Scale = 2.0;
    const double b2Shape = 2.0;
    const double b2Scale = 3.7;
    const double phaseOffset = 4;
    const unsigned DModelParams = 13;
    vec thetaD = { b1Shape, b1Scale, b2Shape, b2Scale, phaseSisterShiftWidths, phaseOffset };


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

    void TestHeModelSPSAFixture()
    {
        double AIC;



        Timer::Reset();
        Timer::Print("Begin TestHeSPSA Fixture @ ");

        //Setup params logfile
        LogFile* p_log = LogFile::Instance();
        p_log->Set(0, "HeModelParameterOptimise", "HeNoise");
        *p_log << "k\tphase2\tphase3\tPP1\tPP2\tPD2\tPP3\n";

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
            unsigned currStartSeed = 0 + 1000*k;
            unsigned currEndSeed = 999 + 1000*k;

            *p_log << k << "\t" << thetaS(0) << "\t" << thetaS(1) << "\t" << thetaS(2) << "\t" << thetaS(3) << "\t"
                    << thetaS(4) << "\t" << thetaS(5) << "\n";

            AIC = evaluateAIC(thetaS, false, DModelParams, numberPoints, currStartSeed, currEndSeed, k);

            *p_log << "AIC: " << AIC << "\n";
            AICvec(k) = AIC;
            k++;

        }

        double AIC_std_dev = stddev(AIC);
        *p_log << "AIC std: " << AIC_std_dev << "\n";

        p_log->Close();
    }

    /*******************
     * UTILITY FUNCTIONS
     *******************/

    double evaluateAIC(vec paramsTheta, bool detMode, unsigned numberParams, unsigned totalPoints, unsigned startSeed, unsigned endSeed, unsigned k)
    {

        //vector to store the 3 residual sums of squares
        vec RSS = vec(3);
        unsigned RSSindex = 0;

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

            uvec histo = hist(runSimulator(inductionTime, paramsTheta, detMode, startSeed, endSeed, k), centers);
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

    uvec runSimulator(double inductionTime, vec theta, bool deterministic, unsigned startSeed, unsigned endSeed, unsigned k)
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
            p_simulator->SetOutputDirectory("UnusedSimOutputOptimisations"); //unused output
            p_simulator->Solve();

            //Count lineage size
            unsigned count = cell_population->GetNumRealCells();
            unsigned seedk = seed - 1000*k;
            counts(seedk) = count;
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
        //if the lineage started before the present induction time, give it a positive TiL and run the full simulation time
        if (lineageStartTime < inductionTime)
        {
            currTiL = inductionTime - lineageStartTime;
            currSimEndTime = endTime - inductionTime;
        }
        //if the lineage starts after the induction time, give it zero TiL run the appropriate-length simulation
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
            p_model->SetModelParameters(currTiL, gammaShift, gammaShape, gammaScale, theta(0), theta(0) + theta(1),
                                        theta(2), phase1PD, theta(3), theta(4), theta(5), phase3PD);
        }
        else
        {
            //Gamma-distribute 8/15hr phase boundaries
            double currPhase2Boundary = theta(5) + p_RNG->GammaRandomDeviate(theta(0), theta(1));
            double currPhase3Boundary = currPhase2Boundary + p_RNG->GammaRandomDeviate(theta(2), theta(3));

            p_model->SetDeterministicMode(currTiL, gammaShift, gammaShape, gammaScale, currPhase2Boundary,
                                          currPhase3Boundary, theta(4));
        }
        return currSimEndTime;
    }

}
;
