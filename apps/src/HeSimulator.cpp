#include <iostream>
#include <string>

#include <cxxtest/TestSuite.h>
#include "ExecutableSupport.hpp"
#include "Exception.hpp"
#include "PetscTools.hpp"
#include "PetscException.hpp"

#include "HeCellCycleModel.hpp"
#include "OffLatticeSimulationPropertyStop.hpp"

#include "AbstractCellBasedTestSuite.hpp"

#include "WildTypeCellMutationState.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

#include "CellsGenerator.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "NodesOnlyMesh.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"

#include "ColumnDataWriter.hpp"

int main(int argc, char *argv[])
{
    ExecutableSupport::StartupWithoutShowingCopyright(&argc, &argv);
    //main() returns code indicating sim run success or failure mode
    int exit_code = ExecutableSupport::EXIT_OK;

    if (argc < 2)
    {
        ExecutableSupport::PrintError(
                "Not enough arguments for simulator. Usage (replace<> with values): HeSimulator 0 <startSeed> <endSeed> <mMitoticModePhase2> <mMitoticModePhase3> <pPP1> <pPD1> <pPP1> <pPD1> <pPP1> <pPD1> <inductionTime> <earliestLineageStart> <latestLineageStart> <endTime> <deterministicModeBool>",
                true);
        exit_code = ExecutableSupport::EXIT_BAD_ARGUMENTS;
    }

    /***********************
     * SIMULATOR PARAMETERS
     ***********************/
    std::string idString, directoryString, filenameString;
    int outputMode; //0 = counts; 1 = mitotic mode events; 2 = mitotic mode sequence sampling
    bool deterministicMode, debugOutput;
    unsigned fixture, startSeed, endSeed; //fixture 0 = He2012; 1 = Wan2016
    double inductionTime, earliestLineageStartTime, latestLineageStartTime, endTime;
    double mitoticModePhase2, mitoticModePhase3, pPP1, pPD1, pPP2, pPD2, pPP3, pPD3; //stochastic model parameters
    double phase1Shape, phase1Scale, phase2Shape, phase2Scale, phaseSisterShiftWidth;

    //PARSE ARGUMENTS
    directoryString = argv[1];
    filenameString = argv[2];
    outputMode = std::stoi(argv[3]);
    deterministicMode = std::stoul(argv[4]);
    fixture = std::stoul(argv[5]);
    debugOutput = std::stoul(argv[6]);
    startSeed = std::stoul(argv[7]);
    endSeed = std::stoul(argv[8]);
    inductionTime = std::stod(argv[9]);
    earliestLineageStartTime = std::stod(argv[10]);
    latestLineageStartTime = std::stod(argv[11]);
    endTime = std::stod(argv[12]);


    if (deterministicMode == 0)
    {
        if (argc < 21)
        {
            ExecutableSupport::PrintError(
                    "Not enough arguments for stochastic model. Usage (replace<> with values): HeSimulator 0 <startSeed> <endSeed> <mMitoticModePhase2> <mMitoticModePhase3> <pPP1> <pPD1> <pPP1> <pPD1> <pPP1> <pPD1> <inductionTime> <earliestLineageStart> <latestLineageStart> <endTime> <deterministicModeBool>",
                    true);
            exit_code = ExecutableSupport::EXIT_BAD_ARGUMENTS;
        }

        mitoticModePhase2 = std::stod(argv[13]);
        mitoticModePhase3 = std::stod(argv[14]);
        pPP1 = std::stod(argv[15]);
        pPD1 = std::stod(argv[16]);
        pPP2 = std::stod(argv[17]);
        pPD2 = std::stod(argv[18]);
        pPP3 = std::stod(argv[19]);
        pPD3 = std::stod(argv[20]);
    }

    else if (deterministicMode == 1)
    {
        if (argc < 18)
        {
            ExecutableSupport::PrintError(
                    "Not enough arguments for deterministic model. Usage (replace<> with values): HeSimulator <startSeed> <endSeed> <mMitoticModePhase2> <mMitoticModePhase3> <pPP1> <pPD1> <pPP1> <pPD1> <pPP1> <pPD1> <inductionTime> <earliestLineageStart> <latestLineageStart> <endTime> <deterministicModeBool>",
                    true);
            exit_code = ExecutableSupport::EXIT_BAD_ARGUMENTS;
        }
        phase1Shape = std::stod(argv[13]);
        phase1Scale = std::stod(argv[14]);
        phase2Shape = std::stod(argv[15]);
        phase2Scale = std::stod(argv[16]);
        phaseSisterShiftWidth = std::stod(argv[17]);
    }

    //Set up singleton LogFile
    LogFile* p_log = LogFile::Instance();
    p_log->Set(0, directoryString, filenameString);

    ExecutableSupport::Print("Simulator writing file " + filenameString + " to directory " + directoryString);

    //Log entry counter
    unsigned entry_number = 1;

    //Write appropriate headers to log
    if (outputMode == 0) *p_log << "Entry\tInduction Time (h)\tSeed\tCount\n";
    if (outputMode == 1) *p_log << "Time (hpf)\tSeed\tCellID\tMitotic Mode (0=PP;1=PD;2=DD)\n";
    if (outputMode == 2) *p_log << "Entry\tSeed\tSequence\n";

    //Instance RNG
    RandomNumberGenerator* p_RNG = RandomNumberGenerator::Instance();

    //Initialise pointers to relevant singleton ProliferativeTypes and Properties
    MAKE_PTR(WildTypeCellMutationState, p_state);
    MAKE_PTR(TransitCellProliferativeType, p_Mitotic);
    MAKE_PTR(DifferentiatedCellProliferativeType, p_PostMitotic);

    //iterate through supplied seed range, executing one simulation per seed
    for (unsigned seed = startSeed; seed <= endSeed; seed++)
    {
        if (outputMode == 2) *p_log << entry_number << "\t" << seed << "\t"; //write seed to log - sequence written by cellcyclemodel objects

        //initialise pointer to debugWriter
        ColumnDataWriter* debugWriter;

        //initialise SimulationTime (permits cellcyclemodel setup)
        SimulationTime::Instance()->SetStartTime(0.0);

        //Reseed the RNG with the required seed
        p_RNG->Reseed(seed);

        //Initialise a HeCellCycleModel and set it up with appropriate TiL values
        HeCellCycleModel* p_cycle_model = new HeCellCycleModel;

        if (debugOutput)
        {
        //Pass ColumnDataWriter to cell cycle model for debug output
            boost::shared_ptr<ColumnDataWriter> p_debugWriter(new ColumnDataWriter(directoryString, filenameString+"DEBUG_"+ std::to_string(seed), false, 10));
            p_cycle_model->EnableModelDebugOutput(p_debugWriter);
            debugWriter = &*p_debugWriter;
        }

        double currTiL; //Time in Lineage offset for lineages induced after first mitosis
        double lineageStartTime; //first mitosis time (hpf)
        double currSimEndTime; //simulation end time (hpf);

        /******************************************************************************
         * Time in Lineage Generation Fixtures & Cell Cycle Model Setup
         ******************************************************************************/

        if (fixture == 0) //He 2012-type fixture - even distribution across nasal-temporal axis
        {
            //generate lineage start time from even random distro across earliest-latest start time figures
            lineageStartTime = (p_RNG->ranf() * (latestLineageStartTime - earliestLineageStartTime))
                    + earliestLineageStartTime;
            //this reflects induction of cells after the lineages' first mitosis
            if (lineageStartTime < inductionTime)
            {
                currTiL = inductionTime - lineageStartTime;
                currSimEndTime = endTime - inductionTime;
                if (outputMode == 1) p_cycle_model->EnableModeEventOutput(inductionTime, seed);
            }
            //if the lineage starts after the induction time, give it zero & TiL run the appropriate-length simulation
            //(ie. the endTime is reduced by the amount of time after induction that the first mitosis occurs)
            if (lineageStartTime >= inductionTime)
            {
                currTiL = 0.0;
                currSimEndTime = endTime - lineageStartTime;
                if (outputMode == 1) p_cycle_model->EnableModeEventOutput(lineageStartTime, seed);
            }

        }
        else if (fixture == 1) //Wan 2016-type fixture - random TiL distribution
        {
            //generate random lineage start time from even random distro across CMZ residency time
            currTiL = p_RNG->ranf() * endTime;
            //calculate simEndTime, only run from tiL to endTime. Min simulation time is 3 minutes (prevents 0 sim endtime errors)
            currSimEndTime = std::max(0.05, endTime - currTiL);
            if (outputMode == 1) p_cycle_model->EnableModeEventOutput(0, seed);
        }

        //Setup lineages' cycle model with appropriate parameters
        p_cycle_model->SetDimension(2);
        p_cycle_model->SetPostMitoticType(p_PostMitotic);

        if (!deterministicMode)
        {
            p_cycle_model->SetModelParameters(currTiL, mitoticModePhase2, mitoticModePhase2 + mitoticModePhase3, pPP1,
                                              pPD1, pPP2, pPD2, pPP3, pPD3);
        }
        else
        {
            //Gamma-distribute phase3 boundary
            double currPhase2Boundary = p_RNG->GammaRandomDeviate(phase1Shape, phase1Scale);
            double currPhase3Boundary = currPhase2Boundary + p_RNG->GammaRandomDeviate(phase2Shape, phase2Scale);

            p_cycle_model->SetDeterministicMode(currTiL, currPhase2Boundary, currPhase3Boundary, phaseSisterShiftWidth);
        }

        //Setup vector containing lineage founder with the properly set up cell cycle model
        std::vector<CellPtr> cells;
        CellPtr p_cell(new Cell(p_state, p_cycle_model));
        p_cell->SetCellProliferativeType(p_Mitotic);
        p_cell->InitialiseCellCycleModel();
        cells.push_back(p_cell);

        //Generate 1x1 mesh for single-cell colony
        HoneycombMeshGenerator generator(1, 1);
        MutableMesh<2, 2>* p_generating_mesh = generator.GetMesh();
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);

        //Setup cell population
        NodeBasedCellPopulation<2>* cell_population(new NodeBasedCellPopulation<2>(mesh, cells));

        //Setup simulator & run simulation
        boost::shared_ptr<OffLatticeSimulationPropertyStop<2>> p_simulator(
                new OffLatticeSimulationPropertyStop<2>(*cell_population));
        p_simulator->SetStopProperty(p_Mitotic); //simulation to stop if no mitotic cells are left
        p_simulator->SetDt(0.05);
        p_simulator->SetEndTime(currSimEndTime);
        p_simulator->SetOutputDirectory("UnusedSimOutput" + filenameString); //unused output
        p_simulator->Solve();

        //Count lineage size
        unsigned count = cell_population->GetNumRealCells();

        if (outputMode == 0) *p_log << entry_number << "\t" << inductionTime << "\t" << seed << "\t" << count << "\n";
        if (outputMode == 2) *p_log << "\n";

        //Reset for next simulation
        SimulationTime::Destroy();
        delete cell_population;
        entry_number++;

        if (debugOutput)
        {
            debugWriter->Close();
        }

    }

    p_RNG->Destroy();
    LogFile::Close();

    return exit_code;
}
;

