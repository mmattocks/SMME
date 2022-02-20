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

    if (argc != 22 && argc != 20)
    {
        ExecutableSupport::PrintError(
                "Wrong arguments for simulator.\nUsage (replace<> with values, pass bools as 0 or 1):\nStochastic Mode:\nHeSimulator <directoryString> <filenameString> <outputModeUnsigned(0=counts,1=events,2=sequence)> <deterministicBool=0> <fixtureUnsigned(0=He;1=Wan;2=test)> <founderAth5Mutant?Bool> <debugOutputBool> <startSeedUnsigned> <endSeedUnsigned>  <inductionTimeDoubleHours> <earliestLineageStartDoubleHours> <latestLineageStartDoubleHours> <endTimeDoubleHours> <mMitoticModePhase2Double> <mMitoticModePhase3Double> <pPP1Double(0-1)> <pPD1Double(0-1)> <pPP1Double(0-1)> <pPD1Double(0-1)> <pPP1Double(0-1)> <pPD1Double(0-1)>\nDeterministic Mode:\nHeSimulator <directoryString> <filenameString> <outputModeUnsigned(0=counts,1=events,2=sequence)> <deterministicBool=1> <fixtureUnsigned(0=He;1=Wan;2=test)> <founderAth5Mutant?Bool> <debugOutputBool> <startSeedUnsigned> <endSeedUnsigned>  <inductionTimeDoubleHours> <earliestLineageStartDoubleHours> <latestLineageStartDoubleHours> <endTimeDoubleHours> <phase1ShapeDouble(>0)> <phase1ScaleDouble(>0)> <phase2ShapeDouble(>0)> <phase2ScaleDouble(>0)> <phaseBoundarySisterShiftWidthDouble>\n",
                true);
        exit_code = ExecutableSupport::EXIT_BAD_ARGUMENTS;
        return exit_code;
    }

    /***********************
     * SIMULATOR PARAMETERS
     ***********************/
    std::string directoryString, filenameString;
    int outputMode; //0 = counts; 1 = mitotic mode events; 2 = mitotic mode sequence sampling
    bool deterministicMode, ath5founder, debugOutput;
    unsigned fixture, startSeed, endSeed; //fixture 0 = He2012; 1 = Wan2016
    double inductionTime, earliestLineageStartTime, latestLineageStartTime, endTime;
    double mitoticModePhase2, mitoticModePhase3, pPP1, pPD1, pPP2, pPD2, pPP3, pPD3; //stochastic model parameters
    double phase1Shape, phase1Scale, phase2Shape, phase2Scale, phaseSisterShiftWidth, phaseOffset;

    //PARSE ARGUMENTS
    directoryString = argv[1];
    filenameString = argv[2];
    outputMode = std::stoi(argv[3]);
    deterministicMode = std::stoul(argv[4]);
    fixture = std::stoul(argv[5]);
    ath5founder = std::stoul(argv[6]);
    debugOutput = std::stoul(argv[7]);
    startSeed = std::stoul(argv[8]);
    endSeed = std::stoul(argv[9]);
    inductionTime = std::stod(argv[10]);
    earliestLineageStartTime = std::stod(argv[11]);
    latestLineageStartTime = std::stod(argv[12]);
    endTime = std::stod(argv[13]);

    if (deterministicMode == 0)
    {
        mitoticModePhase2 = std::stod(argv[14]);
        mitoticModePhase3 = std::stod(argv[15]);
        pPP1 = std::stod(argv[16]);
        pPD1 = std::stod(argv[17]);
        pPP2 = std::stod(argv[18]);
        pPD2 = std::stod(argv[19]);
        pPP3 = std::stod(argv[20]);
        pPD3 = std::stod(argv[21]);
    }
    else if (deterministicMode == 1)
    {
        phase1Shape = std::stod(argv[14]);
        phase1Scale = std::stod(argv[15]);
        phase2Shape = std::stod(argv[16]);
        phase2Scale = std::stod(argv[17]);
        phaseSisterShiftWidth = std::stod(argv[18]);
        phaseOffset = std::stod(argv[19]);
    }
    else
    {
        ExecutableSupport::PrintError("Bad deterministicMode (argument 4). Must be 0 or 1");
        exit_code = ExecutableSupport::EXIT_BAD_ARGUMENTS;
        return exit_code;
    }

    /************************
     * PARAMETER/ARGUMENT SANITY CHECK
     ************************/
    bool sane = 1;

    if (outputMode != 0 && outputMode != 1 && outputMode != 2)
    {
        ExecutableSupport::PrintError(
                "Bad outputMode (argument 3). Must be 0 (counts) 1 (mitotic events) or 2 (sequence sampling)");
        sane = 0;
    }

    if (fixture != 0 && fixture != 1 && fixture != 2)
    {
        ExecutableSupport::PrintError("Bad fixture (argument 5). Must be 0 (He), 1 (Wan), or 2 (validation/test)");
        sane = 0;
    }

    if (ath5founder != 0 && ath5founder != 1)
    {
        ExecutableSupport::PrintError("Bad ath5founder (argument 6). Must be 0 (wild type) or 1 (ath5 mutant)");
        sane = 0;
    }

    if (endSeed < startSeed)
    {
        ExecutableSupport::PrintError("Bad start & end seeds (arguments, 8, 9). endSeed must not be < startSeed");
        sane = 0;
    }

    if (inductionTime >= endTime)
    {
        ExecutableSupport::PrintError("Bad latestLineageStartTime (argument 10). Must be <endTime(arg13)");
        sane = 0;
    }
    if (earliestLineageStartTime >= endTime || earliestLineageStartTime >= latestLineageStartTime)
    {
        ExecutableSupport::PrintError(
                "Bad earliestLineageStartTime (argument 11). Must be <endTime(arg13), <latestLineageStarTime (arg12)");
        sane = 0;
    }
    if (latestLineageStartTime > endTime)
    {
        ExecutableSupport::PrintError("Bad latestLineageStartTime (argument 12). Must be <=endTime(arg13)");
        sane = 0;
    }

    if (deterministicMode == 0)
    {
        if (mitoticModePhase2 < 0)
        {
            ExecutableSupport::PrintError("Bad mitoticModePhase2 (argument 14). Must be >0");
            sane = 0;
        }
        if (mitoticModePhase3 < 0)
        {
            ExecutableSupport::PrintError("Bad mitoticModePhase3 (argument 15). Must be >0");
            sane = 0;
        }
        if (pPP1 + pPD1 > 1 || pPP1 > 1 || pPP1 < 0 || pPD1 > 1 || pPD1 < 0)
        {
            ExecutableSupport::PrintError(
                    "Bad phase 1 probabilities (arguments 16, 17). pPP1 + pPD1 should be >=0, <=1, sum should not exceed 1");
            sane = 0;
        }
        if (pPP2 + pPD2 > 1 || pPP2 > 1 || pPP2 < 0 || pPD2 > 1 || pPD2 < 0)
        {
            ExecutableSupport::PrintError(
                    "Bad phase 2 probabilities (arguments 18, 19). pPP2 + pPD2 should be >=0, <=1, sum should not exceed 1");
            sane = 0;
        }
        if (pPP3 + pPD3 > 1 || pPP3 > 1 || pPP3 < 0 || pPD3 > 1 || pPD3 < 0)
        {
            ExecutableSupport::PrintError(
                    "Bad phase 3 probabilities (arguments 20, 21). pPP3 + pPD3 should be >=0, <=1, sum should not exceed 1");
            sane = 0;
        }
    }

    if (deterministicMode == 1)
    {
        if (phase1Shape <= 0)
        {
            ExecutableSupport::PrintError("Bad phase1Shape (argument 14). Must be >0");
            sane = 0;
        }
        if (phase1Scale <= 0)
        {
            ExecutableSupport::PrintError("Bad phase1Scale (argument 15). Must be >0");
            sane = 0;
        }
        if (phase2Shape <= 0)
        {
            ExecutableSupport::PrintError("Bad phase1Shape (argument 16). Must be >0");
            sane = 0;
        }
        if (phase2Scale <= 0)
        {
            ExecutableSupport::PrintError("Bad phase1Scale (argument 17). Must be >0");
            sane = 0;
        }
        if (phaseSisterShiftWidth <= 0)
        {
            ExecutableSupport::PrintError("Bad phaseSisterShiftWidth (argument 18). Must be >0");
            sane = 0;
        }
    }

    if (sane == 0)
    {
        ExecutableSupport::PrintError("Exiting with bad arguments. See errors for details");
        exit_code = ExecutableSupport::EXIT_BAD_ARGUMENTS;
        return exit_code;

    }

    /************************
     * SIMULATOR OUTPUT SETUP
     ************************/

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
    MAKE_PTR(Ath5Mo, p_Morpholino);
    MAKE_PTR(CellLabel, p_label);

    /************************
     * SIMULATOR SETUP & RUN
     ************************/

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
            boost::shared_ptr<ColumnDataWriter> p_debugWriter(
                    new ColumnDataWriter(directoryString, filenameString + "DEBUG_" + std::to_string(seed), false, 10));
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
            //if the lineage starts after the induction time, give it zero TiL & run the appropriate-length simulation
            //(ie. the endTime is reduced by the amount of time after induction that the first mitosis occurs)
            if (lineageStartTime >= inductionTime)
            {
                currTiL = 0.0;
                currSimEndTime = endTime - lineageStartTime;
                if (outputMode == 1) p_cycle_model->EnableModeEventOutput(lineageStartTime, seed);
            }

        }
        else if (fixture == 1) //Wan 2016-type fixture - each lineage founder selected randomly across residency time, simulator allowed to run until end of residency time
        //passing residency time (as latestLineageStartTime) and endTime separately allows for creation of "shadow CMZ" population
        //this allows investigation of different assumptions about how Wan et al.'s model output was generated
        {
            //generate random lineage start time from even random distro across CMZ residency time
            currTiL = p_RNG->ranf() * latestLineageStartTime;
            currSimEndTime = std::max(.05, endTime - currTiL); //minimum 1 timestep, prevents 0 timestep SimulationTime error
            if (outputMode == 1) p_cycle_model->EnableModeEventOutput(0, seed);
        }
        else if (fixture == 2) //validation fixture- all founders have TiL given by induction time
        {
            currTiL = inductionTime;
            currSimEndTime = endTime;
        }

        //Setup lineages' cycle model with appropriate parameters
        p_cycle_model->SetDimension(2);
        //p_cycle_model->SetPostMitoticType(p_PostMitotic);

        if (!deterministicMode)
        {
            p_cycle_model->SetModelParameters(currTiL, mitoticModePhase2, mitoticModePhase2 + mitoticModePhase3, pPP1,
                                              pPD1, pPP2, pPD2, pPP3, pPD3);
        }
        else
        {
            //Gamma-distribute phase3 boundary
            double currPhase2Boundary = phaseOffset + p_RNG->GammaRandomDeviate(phase1Shape, phase1Scale);
            double currPhase3Boundary = currPhase2Boundary + p_RNG->GammaRandomDeviate(phase2Shape, phase2Scale);

            p_cycle_model->SetDeterministicMode(currTiL, currPhase2Boundary, currPhase3Boundary, phaseSisterShiftWidth);
        }

        if (outputMode == 2) p_cycle_model->EnableSequenceSampler();

        //Setup vector containing lineage founder with the properly set up cell cycle model
        std::vector<CellPtr> cells;
        CellPtr p_cell(new Cell(p_state, p_cycle_model));
        p_cell->SetCellProliferativeType(p_Mitotic);
        if (ath5founder == 1) p_cell->AddCellProperty(p_Morpholino);
        if (outputMode == 2) p_cell->AddCellProperty(p_label);
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
