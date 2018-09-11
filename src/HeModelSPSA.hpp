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

/*************************************************
 * SPSA tuning fixture
 * Intended to estimate appropriate SPSA constants
 *************************************************/

class TestHeModelSPSATune : public AbstractCellBasedTestSuite
{
public:

    /*******************
     * SPSA COEFFICIENTS
     *******************/

    const double a = .45;
    const double c = 1.8; //Approximate standard deviation of 1000-seed AIC measurements at initial thetas
    const double A = 3; //10% of expected # iterations
    const double alpha = .602;
    const double gamma = .101;

    /********************
     * SPSA UTILITY VARS
     ********************/
    const unsigned maxIterations = 30;
    const unsigned numberPoints = 90; //total number of comparison points between model output and empirical data (30 per timepoint)
    const unsigned thetaSize = 4; //4 parameters to optimise
    vec delta = vec(thetaSize); //vector for bernoulli random variables to perturb theta

    /************************
     * SIMULATION PARAMETERS
     ************************/

    //Define start and end RNG seeds; determines:
    //#lineages per loss function run
    //unique sequence of RNG results for each lineage
    const unsigned startSeed = 0;
    const unsigned endSeed = 999;
    unsigned totalSeeds = endSeed - startSeed + 1;

    /************************
     *GLOBAL MODEL PARAMETERS
     ************************/

    //Values defining different marker induction timepoints & relative start time of TiL counter
    const double earliestLineageStartTime = 23.0; //RPCs begin to enter "He model regime" nasally at 23hpf
    const double latestLineageStartTime = 39.0; //last temporal retinal RPC has entered "He model regime" at 39hpf
    const std::vector<double> inductionTimes = { 24.0, 32.0, 48.0 };
    const double endTime = 72.0;
    double currSimEndTime = endTime;
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
    const double mitoticModePhase3 = 15.0;
    const double phase1PP = 1.0;
    const double phase1PD = 0.0;
    const double phase2PP = 0.2;
    const double phase2PD = 0.4;
    const double phase3PP = 0.2;
    const double phase3PD = 0.0;
    const unsigned HeModelParams = 15;
    vec thetaS = { mitoticModePhase2, mitoticModePhase3, phase2PP * 100, phase2PD * 100 };
    vec probScaleVec = { 1, 1, 10, 10 }; //scales ak gain sequence for probability variables
    //DETERMINISTIC MITOTIC MODE
    //Phase boundary shift parameters
    const double phase1Shift = 6.0;
    const double phaseSisterShiftWidths = 0.8;
    const double b1Scale = 1.6;
    const double b2Scale = 3.7;
    const unsigned DModelParams = 11;
    vec thetaD = { phase1Shift, b1Scale, b2Scale, phaseSisterShiftWidths };

    /****************************
     * HE ET AL EMPIRICAL RESULTS
     ****************************/

    const vec he24counts = { 0, 0, 1, 0, 1, 2, 0, 7, 4, 8, 6, 5, 6, 5, 5, 6, 1, 1, 2, 3, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 };
    const vec he32counts = { 7, 20, 25, 21, 24, 17, 11, 7, 8, 5, 7, 2, 1, 6, 5, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0,
                             0, 0 };
    const vec he48counts =
            { 59, 86, 2, 12, 1, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    const unsigned lineagesSampled24 = 64;
    const unsigned lineagesSampled32 = 169;
    const unsigned lineagesSampled48 = 163;

    uvec centers = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27,
                     28, 29, 30 };

    /*********************
     * TEST FIXTURES
     *********************/

    void TestHeModelSPSA()
    {
        Timer::Reset();
        Timer::Print("Begin TestHeSPSA Fixture @ ");

        //Setup params logfile
        LogFile* p_log = LogFile::Instance();
        p_log->Set(0, "HeModelParameterOptimise", "HeParams");
        *p_log << "k\tphase2\tphase3\tPP\tPD\n";

        unsigned k = 0;

        while (k <= maxIterations - 1)
        {
            *p_log << k << "\t" << thetaS(0) << "\t" << thetaS(1) << "\t" << thetaS(2) << "\t" << thetaS(3) << "\n";

            getBernoulliVector(); //populate deltak perturbation vector

            std::stringstream deltareport;
            deltareport << " delta0 " << delta(0) << " delta1 " << delta(1) << " delta2 " << delta(2) << " delta3 "
                    << delta(3);
            Timer::Print(deltareport.str());

            double ak = a / (pow((A + k + 1), alpha)); //calculate ak from gain sequence
            vec scaledak = ak * probScaleVec; // scale ak appropriately for parameters expressed in hrs & percent

            std::stringstream akreport;
            akreport << "ak " << ak << " scaledak0 " << scaledak(0) << " scaledak1 " << scaledak(1) << " scaledak2 "
                    << scaledak(2) << " scaledak 3 " << scaledak(3);
            Timer::Print(akreport.str());

            double ck = c / (pow((k + 1), gamma)); //calculate ck from gain sequence

            vec thetaPlus = thetaS + ck * delta;

            std::stringstream thetaPlusReport;
            thetaPlusReport << "TP0 " << thetaPlus(0) << " TP1 " << thetaPlus(1) << " TP2 " << thetaPlus(2) << " TP3 "
                    << thetaPlus(3);
            Timer::Print(thetaPlusReport.str());

            vec thetaMinus = thetaS - ck * delta;

            std::stringstream thetaMinusReport;
            thetaMinusReport << "TM0 " << thetaMinus(0) << " TM1 " << thetaMinus(1) << " TM2 " << thetaMinus(2)
                    << " TM3 " << thetaMinus(3);
            Timer::Print(thetaMinusReport.str());

            double positiveAIC = evaluateAIC(thetaPlus, false, HeModelParams, numberPoints);
            double negativeAIC = evaluateAIC(thetaMinus, false, HeModelParams, numberPoints);

            std::stringstream AICreport;
            AICreport << " positiveAIC " << positiveAIC << " negativeAIC " << negativeAIC;
            Timer::Print(AICreport.str());

            vec ghat = (positiveAIC - negativeAIC) / 2 * ck * delta;

            std::stringstream ghatreport;
            ghatreport << " ghat0 " << ghat(0) << " ghat1 " << ghat(1) << " ghat2 " << ghat(2) << " ghat 3 " << ghat(3);
            Timer::Print(ghatreport.str());

            std::stringstream oldthetasreport;
            oldthetasreport << " thetas(0) " << thetaS(0) << " thetas(1) " << thetaS(1) << " thetas(2) " << thetaS(2)
                    << " thetas(3) " << thetaS(3);
            Timer::Print(oldthetasreport.str());

            vec iterator = scaledak % ghat;

            std::stringstream iteratorreport;
            iteratorreport << " iterator(0) " << iterator(0) << " iterator(1) " << iterator(1) << " iterator(2) "
                    << iterator(2) << " iterator(3) " << iterator(3);
            Timer::Print(iteratorreport.str());

            thetaS = thetaS - scaledak % ghat;

            //bound off negative values
            if (thetaS(2) < 0) thetaS(2) = 0;
            if (thetaS(3) < 0) thetaS(3) = 0;
            //if, after negative bounding, PP + PD > 1, project appropriately onto legal (PP,PD) parameter space
            if (thetaS(2) + thetaS(3) > 1)
            {
                if (thetaS(2) > 1 && thetaS(3) < thetaS(2) - 1) // these (PP,PD) points will project into negative PD space
                {
                    thetaS(2) = 1;
                    thetaS(3) = 0;
                }
                if (thetaS(3) > 1 && thetaS(2) < thetaS(3) - 1) // these (PP,PD) points will project into negative PP space
                {
                    thetaS(3) = 1;
                    thetaS(2) = 0;
                }

                /******************************
                 * SPSA PROBABILITY THETA CONSTRAINTS
                 ******************************/

                std::valarray<double> p1 = { thetaS(2), thetaS(3) }; // compose (PP,PD) valarray point
                std::valarray<double> v1 = { 1, -1 }; //vector from (0,1) to (1,0)
                std::valarray<double> v2 = { thetaS(2), thetaS(3) - 1 }; //vector from (0,1) to p1
                double dotProduct = (v1 * v2).sum();
                double lengthv1 = sqrt(2);
                double lengthv2 = sqrt(pow(thetaS(2), 2) + pow(thetaS(3) - 1, 2));
                double cos = dotProduct / (lengthv1 * lengthv2);
                double projLenOfLine = cos * lengthv2;
                double newPP = projLenOfLine / lengthv1;
                double newPD = 1 - projLenOfLine / lengthv2;

                std::stringstream newthetasreport;
                newthetasreport << " thetas(0) " << thetaS(0) << " thetas(1) " << thetaS(1) << " thetas(2) "
                        << thetaS(2) << " thetas(3) " << thetaS(3);
                Timer::Print(newthetasreport.str());

                k++;
            }
            p_log->Close();
        }
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

    double evaluateAIC(vec paramsTheta, bool detMode, unsigned numberParams, unsigned totalPoints)
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
                empiricalCounts = he24counts;
                empiricalLineagesSampled = lineagesSampled24;
            }
            if (inductionTime == 32)
            {
                empiricalCounts = he32counts;
                empiricalLineagesSampled = lineagesSampled32;
            }
            if (inductionTime == 48)
            {
                empiricalCounts = he48counts;
                empiricalLineagesSampled = lineagesSampled48;
            }

            uvec histo = hist(runSimulator(inductionTime, paramsTheta, detMode), centers);
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
            p_simulator->SetOutputDirectory("UnusedSimOutputOptimisations"); //unused output
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
            p_model->SetModelParameters(currTiL, gammaShift, gammaShape, gammaScale, theta(0), theta(1), phase1PP,
                                        phase1PD, theta(2) / 100, theta(3) / 100, phase3PP, phase3PD);
        }
        else
        {
            //Gamma-distribute 8/15hr phase boundaries
            double currPhase2Boundary = theta(0) + p_RNG->GammaRandomDeviate(gammaShape, theta(1));
            double currPhase3Boundary = currPhase2Boundary + p_RNG->GammaRandomDeviate(gammaShape, theta(2));

            p_model->SetDeterministicMode(currTiL, gammaShift, gammaShape, gammaScale, currPhase2Boundary,
                                          currPhase3Boundary, theta(3));
        }
        return currSimEndTime;
    }

}
;
