#include <cxxtest/TestSuite.h>
#include "Exception.hpp"
#include <string>

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



class TestHeModelInduction48 : public AbstractCellBasedTestSuite
{
public:
    bool debugOutput = false;

    /************************
     * SIMULATION PARAMETERS
     ************************/

    //Define start and end RNG seeds; determines:
    //#lineages
    //unique sequence of RNG results for each lineage
    const unsigned startSeed = 20000;
    const unsigned endSeed = 22999;

    /************************
     *GLOBAL MODEL PARAMETERS
     ************************/

    //Values defining different marker induction timepoints & relative start time of TiL counter
    const double earliestLineageStartTime = 23.0; //RPCs begin to enter "He model regime" nasally at 23hpf
    const double latestLineageStartTime = 39.0; //last temporal retinal RPC has entered "He model regime" at 39hpf
    const double inductionTime = 48.0;
    const std::string inductionTimeString = "48";
    const double endTime = 72.0;
    double currSimEndTime = endTime;
    //Cell cycle length gamma PDF parameters
    const double gammaShift = 4.0;
    const double gammaShape = 2.0;
    const double gammaScale = 1.0;
    //3-phase mitotic mode time periodisation

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

    void TestHeInductionCount48Fixture()
    {
        Timer::Reset();
        Timer::Print("Begin TestHeInductionCount" + inductionTimeString + "Fixture @ ");

        //Initialise pointers to relevant ProliferativeTypes and Properties
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_Mitotic);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_PostMitotic);

        //Setup count logfile
        LogFile* p_log = LogFile::Instance();
        p_log->Set(0, "HeInductionCounts", "He" + inductionTimeString);
        *p_log << "Entry\tInduction Time (h)\tSeed\tCount\n";

        //Setup RNG
        RandomNumberGenerator* p_RNG = RandomNumberGenerator::Instance();

        //iterate through supplied seed range, executing one simulation per seed
        for (unsigned seed = startSeed; seed <= endSeed; seed++)
        {
            //initialise pointer to debugWriter
            ColumnDataWriter* debugWriter;

            //entry counter for count writer
            unsigned entry_number = 1;

            //Reseed the RNG with the required seed
            p_RNG->Reseed(seed);

            //Initialise a HeCellCycleModel and set it up with appropriate TiL values, also obtaining the correct simulation end time
            HeCellCycleModel* p_cycle_model = new HeCellCycleModel;
            double simEndTime = SetupNTHeModel(p_cycle_model, p_PostMitotic);

            if (debugOutput)
            {
                //Setup the writer and model for debug output
                debugWriter = SetupDebugOutput(p_cycle_model, inductionTime, "He", seed);
            }

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
            p_simulator->SetOutputDirectory("UnusedSimOutputHe" + inductionTimeString); //unused output
            p_simulator->Solve();

            //Count lineage size
            unsigned count = cell_population->GetNumRealCells();
            double currSimTime = SimulationTime::Instance()->GetTime();

            //Output simulation report to console
            std::stringstream reportStream;
            reportStream << "Simulation Seed " << seed << ": Final Lineage Count: " << count
                    << " Sim Complete; Last SimTime " << currSimTime << " @ ";
            auto simReport = reportStream.str();
            Timer::Print(simReport);

            //Write simulation report to logfile
            *p_log << entry_number << "\t" << inductionTime << "\t" << seed << "\t" << count << "\n";

            //Reset for next simulation
            entry_number++;
            SimulationTime::Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);
            delete cell_population;

            if (debugOutput)
            {
                debugWriter->Close();
            }
        }
        p_RNG->Destroy();
        LogFile::Close();
    }

    void TestDeterministicModelInduction48Fixture()
    {
        Timer::Reset();
        Timer::Print("Begin TestDeterministicModelInduction" + inductionTimeString + "Fixture @ ");

        //Initialise pointers to relevant ProliferativeTypes and Properties
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_Mitotic);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_PostMitotic);

        //Setup count logfile
        LogFile* p_log = LogFile::Instance();
        p_log->Set(0, "HeInductionCounts", "Deterministic" + inductionTimeString);
        *p_log << "Entry\tInduction Time (h)\tSeed\tCount\n";

        //Setup RNG
        RandomNumberGenerator* p_RNG = RandomNumberGenerator::Instance();

        //iterate through supplied seed range, executing one simulation per seed
        for (unsigned seed = startSeed; seed <= endSeed; seed++)
        {
            //initialise pointer to debugWriter
            ColumnDataWriter* debugWriter;

            //entry counter for count writer
            unsigned entry_number = 1;

            //Reseed the RNG with the required seed
            p_RNG->Reseed(seed);

            //Initialise a deterministic HeCellCycleModel and set it up with appropriate TiL values, also obtaining the correct simulation end time
            HeCellCycleModel* p_cycle_model = new HeCellCycleModel;
            double simEndTime = SetupNTHeModel(p_cycle_model, p_PostMitotic, true);

            if (debugOutput)
            {
                //Setup the writer and model for debug output
                debugWriter = SetupDebugOutput(p_cycle_model, inductionTime, "Deterministic", seed);
            }

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
            p_simulator->SetOutputDirectory("UnusedSimOutputDeterministicInduction" + inductionTimeString);
            p_simulator->Solve();

            //Count lineage size
            unsigned count = cell_population->GetNumRealCells();
            double currSimTime = SimulationTime::Instance()->GetTime();

            //Output simulation report to console
            std::stringstream reportStream;
            reportStream << "Simulation Seed " << seed << ": Final Lineage Count: " << count
                    << " Sim Complete; Last SimTime " << currSimTime << " @ ";
            auto simReport = reportStream.str();
            Timer::Print(simReport);

            //Write simulation report to logfile
            *p_log << entry_number << "\t" << inductionTime << "\t" << seed << "\t" << count << "\n";

            //Reset for next simulation
            entry_number++;
            SimulationTime::Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);
            delete cell_population;

            if (debugOutput)
            {
                debugWriter->Close();
            }
        }
        p_RNG->Destroy();
        LogFile::Close();
    }

    /*******************
     * UTILITY FUNCTIONS
     *******************/

    double SetupNTHeModel(HeCellCycleModel* p_model, boost::shared_ptr<AbstractCellProperty> pmType,
                          bool deterministic = false, std::string debugDir = "debug")
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
