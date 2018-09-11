#include <cxxtest/TestSuite.h>
#include "Exception.hpp"
#include <string>
#include <iomanip>

#include "../../projects/ISP/src/HeCellCycleModel.hpp"
#include "../../projects/ISP/src/OffLatticeSimulationPropertyStop.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "PetscSetupAndFinalize.hpp"

#include "AbstractCellBasedTestSuite.hpp"

#include "HoneycombMeshGenerator.hpp"
#include "NodesOnlyMesh.hpp"

#include "CellsGenerator.hpp"
#include "WildTypeCellMutationState.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"

#include "SmartPointers.hpp"

#include "ColumnDataWriter.hpp"

class TestHeModelModeRate80 : public AbstractCellBasedTestSuite
{
public:
    bool debugOutput = false;

    /************************
     * SIMULATION PARAMETERS
     ************************/

    //Define start and end RNG seeds; determines:
    //#lineages
    //reproducible sequence of RNG results for each lineage
    const unsigned startSeed = 0;
    const unsigned endSeed = 2500;

    /************************
     *GLOBAL MODEL PARAMETERS
     ************************/

    //Values defining different marker induction timepoints & relative start time of TiL counter
    const double earliestLineageStartTime = 23.0; //RPCs begin to enter "He model regime" nasally at 23hpf
    const double latestLineageStartTime = 39.0; //last temporal retinal RPC has entered "He model regime" at 39hpf
    const double inductionTime = 23.0;
    const std::string inductionTimeString = "23";
    const double endTime = 80.0;
    double currSimEndTime = endTime;
    //Cell cycle length gamma PDF parameters
    const double gammaShift = 4.0;
    const double gammaShape = 2.0;
    const double gammaScale = 1.0;

    /**************************
     * SPECIFIC MODEL PARAMETERS
     **************************/
    //STOCHASTIC MITOTIC MODE
    //3-phase mitotic mode time periodisation
    const double mitoticModePhase2 = 2.798;
    const double mitoticModePhase3 = 18.307;
    //mitotic mode per-phase probabilities
    const double phase1PP = 1.0;
    const double phase1PD = 0.0;
    const double phase2PP = 0.3685;
    const double phase2PD = 0.4811;
    const double phase3PP = 0.2;
    const double phase3PD = 0.0;
    //DETERMINISTIC MITOTIC MODE
    //Phase boundary shift parameters
    const double phase1Shift = .383;
    const double phaseSisterShiftWidths = 0.789;
    const double b1Scale = 1.879;
    const double b2Scale = 9.226;
    /************************
     * PARAMETER SANITY CHECK
     ************************/
    void TestParameters()
    {
        TS_ASSERT_LESS_THAN(mitoticModePhase2, mitoticModePhase3);
        TS_ASSERT_LESS_THAN(earliestLineageStartTime, endTime);
        TS_ASSERT_LESS_THAN(latestLineageStartTime, endTime);
        TS_ASSERT_LESS_THAN(inductionTime, endTime);
        TS_ASSERT_LESS_THAN(earliestLineageStartTime, latestLineageStartTime);
        TS_ASSERT_EQUALS(1.0, phase1PP + phase1PD + (1.0 - phase1PP - phase1PD));
        TS_ASSERT_EQUALS(1.0, phase2PP + phase2PD + (1.0 - phase2PP - phase2PD));
        TS_ASSERT_EQUALS(1.0, phase3PP + phase3PD + (1.0 - phase3PP - phase3PD));
    }

    /*********************
     * TEST FIXTURES
     *********************/

    void TestHeMitoticModeRateFixture()
    {
        Timer::Reset();
        Timer::Print("Begin TestHeMitoticModeRateFixture @ ");

        //Initialise pointers to relevant ProliferativeTypes and Properties
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_Mitotic);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_PostMitotic);

        //Setup event logfile
        LogFile* p_log = LogFile::Instance();
        p_log->Set(0, "HeMitoticRateEvents", "He" + inductionTimeString);
        *p_log << "Time (hpf)\tSeed\tCellID\tMitotic Mode (0=PP;1=PD;2=DD)\n";

        //Setup RNG
        RandomNumberGenerator* p_RNG = RandomNumberGenerator::Instance();

        //iterate through supplied seed range, executing one simulation per seed
        for (unsigned seed = startSeed; seed <= endSeed; seed++)
        {
            //Reseed the RNG with the required seed
            p_RNG->Reseed(seed);

            //Setup lineages' cycle model with appropriate parameters
            HeCellCycleModel* p_cycle_model = new HeCellCycleModel();
            double simEndTime = SetupNTHeModelEvents(p_cycle_model, p_PostMitotic, seed);

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
            p_simulator->SetEndTime(simEndTime); //end simulations at 80hr TiL
            p_simulator->SetOutputDirectory("UnusedSimOutputHeMitoticModeRateFixture");
            p_simulator->Solve();

            unsigned count = cell_population->GetNumRealCells();
            double currSimTime = SimulationTime::Instance()->GetTime();

            //Output simulation report
            std::stringstream reportStream;
            reportStream << "Simulation Seed " << seed << ": Final Lineage Count: " << count
                    << " Sim Complete; Last SimTime " << currSimTime << " @ ";
            auto simReport = reportStream.str();
            Timer::Print(simReport);

            //Reset for next simulation
            SimulationTime::Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);
            delete cell_population;
        }

        p_RNG->Destroy();
        LogFile::Close();
    }

    void TestDeterministicModelModeRateFixture()
    {
        Timer::Reset();
        Timer::Print("Begin TestDeterministicModelModeRateFixture @ ");

        //Initialise pointers to relevant ProliferativeTypes and Properties
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_Mitotic);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_PostMitotic);

        //Setup event logfile
        LogFile* p_log = LogFile::Instance();
        p_log->Set(0, "HeMitoticRateEvents", "Deterministic" + inductionTimeString);
        *p_log << "Entry\tInduction Time (h)\tSeed\tCount\n";

        //Setup RNG
        RandomNumberGenerator* p_RNG = RandomNumberGenerator::Instance();

        //iterate through supplied seed range, executing one simulation per seed
        for (unsigned seed = startSeed; seed <= endSeed; seed++)
        {
            //Reseed the RNG with the required seed
            p_RNG->Reseed(seed);

            //Setup lineages' cycle model with appropriate parameters
            HeCellCycleModel* p_cycle_model = new HeCellCycleModel();
            double simEndTime = SetupNTHeModelEvents(p_cycle_model, p_PostMitotic, seed, true);

            //Generate 1x1 mesh for single-cell colony
            HoneycombMeshGenerator generator(1, 1);
            MutableMesh<2, 2>* p_generating_mesh = generator.GetMesh();
            //Prepare nodes-only mesh
            NodesOnlyMesh<2> mesh;
            mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);

            //Setup vector containing lineage founder*/
            std::vector<CellPtr> cells;
            CellPtr p_cell(new Cell(p_state, p_cycle_model));
            p_cell->SetCellProliferativeType(p_Mitotic);
            p_cell->InitialiseCellCycleModel();
            cells.push_back(p_cell);

            //Setup cell population
            NodeBasedCellPopulation<2>* cell_population(new NodeBasedCellPopulation<2>(mesh, cells));

            //Setup simulator & run simulation

            boost::shared_ptr<OffLatticeSimulationPropertyStop<2>> p_simulator(
                    new OffLatticeSimulationPropertyStop<2>(*cell_population));
            p_simulator->SetStopProperty(p_Mitotic); //simulation to stop if no mitotic cells are left
            p_simulator->SetDt(0.05);
            p_simulator->SetEndTime(simEndTime);
            p_simulator->SetOutputDirectory("UnusedSimOutputDeterministicMitoticModeRateFixture");
            p_simulator->Solve();

            unsigned count = cell_population->GetNumRealCells();
            double currSimTime = SimulationTime::Instance()->GetTime();

            //Output simulation report
            std::stringstream reportStream;
            reportStream << "Simulation Seed " << seed << ": Final Lineage Count: " << count
                    << " Sim Complete; Last SimTime " << currSimTime << " @ ";
            auto simReport = reportStream.str();
            Timer::Print(simReport);

            //Reset for next simulation
            SimulationTime::Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);
            delete cell_population;

        }

        p_RNG->Destroy();
        LogFile::Close();
    }

    /*******************
     * UTILITY FUNCTIONS
     *******************/

    double SetupNTHeModelEvents(HeCellCycleModel* p_model, boost::shared_ptr<AbstractCellProperty> pmType,
                                unsigned seed, bool deterministic = false, std::string debugDir = "debug")
    {
        //get the RNG instance
        RandomNumberGenerator* p_RNG = RandomNumberGenerator::Instance();
        /******************************************************
         * Generate Time in Lineage figure from even distribution across nasal-temporal gradient
         ******************************************************/
        //if the inductionTime is less than the latestLineageStartTime
        //generate random lineage start time from even random distro across earliest-latest start time figures
        double lineageStartTime = (p_RNG->ranf() * (latestLineageStartTime - earliestLineageStartTime))
                + earliestLineageStartTime;
        double eventWriterStartTime;
        double currTiL;
        double currSimEndTime;
        //if the lineage started before the present induction time, give it a positive TiL and run the full simulation time
        if (lineageStartTime < inductionTime)
        {
            currTiL = inductionTime - lineageStartTime;
            currSimEndTime = endTime - inductionTime;
            eventWriterStartTime = inductionTime;
        }
        //if the lineage starts after the induction time, give it zero TiL run the appropriate-length simulation
        if (lineageStartTime >= inductionTime)
        {
            currTiL = 0.0;
            currSimEndTime = endTime - lineageStartTime;
            eventWriterStartTime = lineageStartTime;
        }

        //Setup lineages' cycle model with appropriate parameters
        p_model->SetDimension(2);
        p_model->SetPostMitoticType(pmType);
        p_model->EnableModeEventOutput(eventWriterStartTime, seed);

        if (!deterministic)
        {
            p_model->SetModelParameters(currTiL, gammaShift, gammaShape, gammaScale, mitoticModePhase2,
                                        mitoticModePhase3, phase1PP, phase1PD, phase2PP, phase2PD, phase3PP, phase3PD);
        }
        else
        {
            //Gamma-distribute 8/15hr phase boundaries
            double currPhase2Boundary = phase1Shift + p_RNG->GammaRandomDeviate(gammaShape, b1Scale);
            double currPhase3Boundary = currPhase2Boundary + p_RNG->GammaRandomDeviate(gammaShape, b2Scale);

            p_model->SetDeterministicMode(currTiL, gammaShift, gammaShape, gammaScale, currPhase2Boundary,
                                          currPhase3Boundary, phaseSisterShiftWidths);
        }
        return currSimEndTime;
    }

    ColumnDataWriter* SetupDebugOutput(HeCellCycleModel* p_model, double time, std::string model, unsigned seed)
    {
        //Pass ColumnDataWriter to cell cycle model for debug output
        std::string debugDir = "TestHeModelInductionDebugFiles";
        std::stringstream filenameStream;
        filenameStream << "HeModelDebug" << model << "_I_" << time << "_S_" << seed;
        std::string debugFilename = filenameStream.str();
        boost::shared_ptr<ColumnDataWriter> p_debugWriter(new ColumnDataWriter(debugDir, debugFilename, false, 10));
        p_model->EnableModelDebugOutput(p_debugWriter);
        ColumnDataWriter* debugWriter = &*p_debugWriter;
        return debugWriter;
    }
};

