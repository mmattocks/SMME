#include <cxxtest/TestSuite.h>
#include "Exception.hpp"
#include <string>
#include <iomanip>

#include "../../projects/ISP/src/HeAth5Mo.hpp"
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

class TestHeModelAth5 : public AbstractCellBasedTestSuite
{
public:

    bool debugOutput = false;

    /************************
     * SIMULATION PARAMETERS
     ************************/

    //Define start and end RNG seeds; determines:
    //#lineages
    //unique sequence of RNG results for each lineage
    const unsigned startSeed = 0;
    const unsigned endSeed = 2500;

    /************************
     *GLOBAL MODEL PARAMETERS
     ************************/

    //Values defining different marker induction timepoints & relative start time of TiL counter
    const double inductionTime = 0.0;
    const std::string inductionTimeString = "0";
    const double endTime = 120.0;
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
        TS_ASSERT_LESS_THAN(inductionTime, endTime);
        TS_ASSERT_EQUALS(1.0, phase1PP + phase1PD + (1.0 - phase1PP - phase1PD));
        TS_ASSERT_EQUALS(1.0, phase2PP + phase2PD + (1.0 - phase2PP - phase2PD));
        TS_ASSERT_EQUALS(1.0, phase3PP + phase3PD + (1.0 - phase3PP - phase3PD));
    }

    /*********************
     * TEST FIXTURES
     *********************/

    void TestHeModelAth5Fixture()
    {
        Timer::Reset();
        Timer::Print("Begin TestHeModelAth5Fixture @ ");

        //Initialise pointers to relevant ProliferativeTypes and Properties
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_Mitotic);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_PostMitotic);
        MAKE_PTR(Ath5Mo, p_Morpholino);

        //Count Writer Setup
        ColumnDataWriter* countWriter = SetupCountWriter("TestHeInductionCounts", "Ath5HeModelResults");

        //Event Writer Initialisation
        boost::shared_ptr<ColumnDataWriter> p_eventWriter (new ColumnDataWriter("HeModelMitoticRateEvents", "Ath5HeModelEvents", false, 10));
        std::vector<int> eventWriterVarIDs = SetupEventWriter(p_eventWriter);

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
            SetupValidateHeModelEvents(p_cycle_model, p_PostMitotic, seed, p_eventWriter, eventWriterVarIDs);

            if (debugOutput)
            {
                //Setup the writer and model for debug output
                debugWriter = SetupDebugOutput(p_cycle_model, inductionTime, "HeAth5", seed);
            }

            //Setup vector containing lineage founder*/
            std::vector<CellPtr> cells;
            CellPtr p_cell(new Cell(p_state, p_cycle_model));
            p_cell->SetCellProliferativeType(p_Mitotic);
            p_cell->AddCellProperty(p_Morpholino);
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
            p_simulator->SetEndTime(endTime);
            p_simulator->SetOutputDirectory("UnusedSimOutputAth5HeModelFixture");
            p_simulator->Solve();

            //Count lineage size
            unsigned count = cell_population->GetNumRealCells();
            double currSimTime = SimulationTime::Instance()->GetTime();

            //Output simulation report
            std::stringstream reportStream;
            reportStream << "Simulation Seed " << seed << ": Final Lineage Count: " << count
                    << " Sim Complete; Last SimTime " << currSimTime << " @ ";
            auto simReport = reportStream.str();
            Timer::Print(simReport);

            //Write simulation report to file
            WriteCountReport(countWriter, entry_number, inductionTime, seed, count);

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
        p_eventWriter->Close();
        countWriter->Close();
    }

    void TestDeterministicModelAth5Fixture()
    {
        Timer::Reset();
        Timer::Print("Begin TestDeterministicModelAth5Fixture @ ");

        //Initialise pointers to relevant ProliferativeTypes and Properties
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_Mitotic);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_PostMitotic);
        MAKE_PTR(Ath5Mo, p_Morpholino);

        //Count Writer Setup
        ColumnDataWriter* countWriter = SetupCountWriter("TestHeInductionCounts", "Ath5DeterministicModelResults");

        //Event Writer Initialisation
        boost::shared_ptr<ColumnDataWriter> p_eventWriter (new ColumnDataWriter("HeModelMitoticRateEvents", "Ath5DeterministicModelMitoticRateEvents", false, 10));
        std::vector<int> eventWriterVarIDs = SetupEventWriter(p_eventWriter);

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
            SetupValidateHeModelEvents(p_cycle_model, p_PostMitotic, seed, p_eventWriter, eventWriterVarIDs, true);

            if (debugOutput)
            {
                //Setup the writer and model for debug output
                debugWriter = SetupDebugOutput(p_cycle_model, inductionTime, "DeterministicAth5", seed);
            }

            //Setup vector containing lineage founder*/
            std::vector<CellPtr> cells;
            CellPtr p_cell(new Cell(p_state, p_cycle_model));
            p_cell->SetCellProliferativeType(p_Mitotic);
            p_cell->AddCellProperty(p_Morpholino);
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
            p_simulator->SetEndTime(endTime);
            p_simulator->SetOutputDirectory("UnusedSimOutputAth5DeterministicModelFixture");
            p_simulator->Solve();

            //Count lineage size
            unsigned count = cell_population->GetNumRealCells();
            double currSimTime = SimulationTime::Instance()->GetTime();

            //Output simulation report
            std::stringstream reportStream;
            reportStream << "Simulation Seed " << seed << ": Final Lineage Count: " << count
                    << " Sim Complete; Last SimTime " << currSimTime << " @ ";
            auto simReport = reportStream.str();
            Timer::Print(simReport);

            //Write simulation report to file
            WriteCountReport(countWriter, entry_number, inductionTime, seed, count);

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
        p_eventWriter->Close();
        countWriter->Close();
    }

    void TestShiftedDeterministicModelAth5Fixture()
    {
        Timer::Reset();
        Timer::Print("Begin TestDeterministicModelAth5Fixture @ ");

        //Initialise pointers to relevant ProliferativeTypes and Properties
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_Mitotic);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_PostMitotic);
        MAKE_PTR(Ath5Mo, p_Morpholino);
        
        //Count Writer Setup
        ColumnDataWriter* countWriter = SetupCountWriter("TestHeInductionCounts", "Ath5ShiftedDeterministicModelResults");

        //Event Writer Initialisation
        boost::shared_ptr<ColumnDataWriter> p_eventWriter (new ColumnDataWriter("HeModelMitoticRateEvents", "Ath5ShiftedDeterministicMitoticRateEvents", false, 10));
        std::vector<int> eventWriterVarIDs = SetupEventWriter(p_eventWriter);

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

            //Initialise a shifted deterministic HeCellCycleModel and set it up with appropriate TiL values, also obtaining the correct simulation end time
            HeCellCycleModel* p_cycle_model = new HeCellCycleModel;
            SetupValidateHeModelEvents(p_cycle_model, p_PostMitotic, seed, p_eventWriter, eventWriterVarIDs, true, true);

            if (debugOutput)
            {
                //Setup the writer and model for debug output
                debugWriter = SetupDebugOutput(p_cycle_model, inductionTime, "ShiftedDeterministicAth5", seed);
            }

            //Setup vector containing lineage founder*/
            std::vector<CellPtr> cells;
            CellPtr p_cell(new Cell(p_state, p_cycle_model));
            p_cell->SetCellProliferativeType(p_Mitotic);
            p_cell->AddCellProperty(p_Morpholino);
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
            p_simulator->SetEndTime(endTime);
            p_simulator->SetOutputDirectory("UnusedSimOutputShiftedAth5DeterministicModelFixture");
            p_simulator->Solve();

            //Count lineage size
            unsigned count = cell_population->GetNumRealCells();
            double currSimTime = SimulationTime::Instance()->GetTime();

            //Output simulation report
            std::stringstream reportStream;
            reportStream << "Simulation Seed " << seed << ": Final Lineage Count: " << count
                    << " Sim Complete; Last SimTime " << currSimTime << " @ ";
            auto simReport = reportStream.str();
            Timer::Print(simReport);

            //Write simulation report to file
            WriteCountReport(countWriter, entry_number, inductionTime, seed, count);

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
        p_eventWriter->Close();
        countWriter->Close();
    }

    /*******************
     * UTILITY FUNCTIONS
     *******************/

    ColumnDataWriter* SetupCountWriter(std::string writeDirectory, std::string filename)
    {
        /************************
         *COUNT WRITER SETUP
         ************************/
        const std::string time = "InductionTime";
        const std::string timeUnit = "h";
        const std::string entry = "Entry";
        const std::string seed = "Seed";
        const std::string count = "Count";
        const std::string number = "No";

        ColumnDataWriter* dataWriter = new ColumnDataWriter(writeDirectory, filename, false, 10);
        dataWriter->DefineUnlimitedDimension(entry, number);
        dataWriter->DefineVariable(time, timeUnit);
        dataWriter->DefineVariable(seed, number);
        dataWriter->DefineVariable(count, number);

        dataWriter->EndDefineMode();
        return (dataWriter);
    }

    void WriteCountReport(ColumnDataWriter* writer, unsigned number, double time, unsigned seed, unsigned count)
    {
        writer->PutVariable(-1, number);
        writer->PutVariable(0, time);
        writer->PutVariable(1, seed);
        writer->PutVariable(2, count);
        writer->AdvanceAlongUnlimitedDimension();
    }

    ColumnDataWriter* SetupDebugOutput(HeCellCycleModel* p_model, double time, std::string model, unsigned seed)
    {
        //Pass ColumnDataWriter to cell cycle model for debug output
        std::string debugDir = "TestHeModelValidationDebugFiles";
        std::stringstream filenameStream;
        filenameStream << "HeModelDebug" << model << "_I_" << time << "_S_" << seed;
        std::string debugFilename = filenameStream.str();
        boost::shared_ptr<ColumnDataWriter> p_debugWriter(
                new ColumnDataWriter(debugDir, debugFilename, false, 10));
        p_model->EnableModelDebugOutput(p_debugWriter);
        ColumnDataWriter* debugWriter = &*p_debugWriter;
        return debugWriter;
    }
    
        std::vector<int> SetupEventWriter(boost::shared_ptr<ColumnDataWriter> writer)
    {
        /*********************
         * EVENT WRITER SETUP
         *********************/
        std::string time = "Time";
        std::string timeUnit = "h";
        std::string seedString = "Seed";
        std::string cellID = "CellID";
        std::string number = "No";
        std::string mitoticMode = "MitoticMode"; //0=PP;1=PD;2=DD
        std::string mode = "Mode";

        std::vector<int> varIDs;

        writer->DefineUnlimitedDimension(time, timeUnit);
        varIDs.push_back(writer->DefineVariable(seedString, number));
        varIDs.push_back(writer->DefineVariable(cellID, number));
        varIDs.push_back(writer->DefineVariable(mitoticMode, mode));

        //Mode comment
        std::string modeComment = "Mitotic Mode: 0 = PP; 1 = PD; 2 = DD";
        writer->SetCommentForInfoFile(modeComment);

        writer->EndDefineMode();

        return varIDs;
    }
    
    void SetupValidateHeModelEvents(HeCellCycleModel* p_model, boost::shared_ptr<AbstractCellProperty> pmType, unsigned seed, boost::shared_ptr<ColumnDataWriter> writer, std::vector<int> varIDs, bool deterministic = false, bool shift = false, std::string debugDir = "debug")
    {
        //get the RNG instance
        RandomNumberGenerator* p_RNG = RandomNumberGenerator::Instance();

        //Setup lineages' cycle model with appropriate parameters
        p_model->SetDimension(2);
        p_model->SetPostMitoticType(pmType);
        p_model->EnableModeEventOutput(inductionTime, seed);

        if (!deterministic)
        {
            p_model->SetModelParameters(inductionTime, gammaShift, gammaShape, gammaScale, mitoticModePhase2,
                                        mitoticModePhase3, phase1PP, phase1PD, phase2PP, phase2PD, phase3PP, phase3PD);
        }
        else
        {
            if (shift)
            {
                //Gamma-distribute 8/15hr phase boundaries
                double currPhase2Boundary = phase1Shift + p_RNG->GammaRandomDeviate(gammaShape, b1Scale);
                double currPhase3Boundary = currPhase2Boundary + p_RNG->GammaRandomDeviate(gammaShape, b2Scale);

                p_model->SetDeterministicMode(inductionTime, gammaShift, gammaShape, gammaScale, currPhase2Boundary,
                                              currPhase3Boundary, phaseSisterShiftWidths);
            }
            else
            {
                p_model->SetDeterministicMode(inductionTime, gammaShift, gammaShape, gammaScale, mitoticModePhase2,
                                              mitoticModePhase3, phaseSisterShiftWidths);
            }
        }
    }
};
