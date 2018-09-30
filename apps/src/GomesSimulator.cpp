#include <iostream>
#include <string>

#include <cxxtest/TestSuite.h>
#include "ExecutableSupport.hpp"
#include "Exception.hpp"
#include "PetscTools.hpp"
#include "PetscException.hpp"

#include "GomesCellCycleModel.hpp"
#include "OffLatticeSimulationPropertyStop.hpp"

#include "AbstractCellBasedTestSuite.hpp"

#include "WildTypeCellMutationState.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "GomesRetinalNeuralFates.hpp"

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

    if (argc != 15)
    {
        ExecutableSupport::PrintError(
                "Wrong arguments for simulator.\nUsage (replace<> with values, pass bools as 0 or 1):\n GomesSimulator <directoryString> <filenameString> <outputModeUnsigned(0=counts,1=events,2=sequence)> <debugOutputBool> <startSeedUnsigned> <endSeedUnsigned> <endTimeDoubleHours> <cellCycleNormalMeanDouble> <cellCycleNormalStdDouble> <pPPDouble(0-1)> <pPDDouble(0-1)> <pBCDouble(0-1)> <pACDouble(0-1)> <pMGDouble(0-1)>",
                true);
        exit_code = ExecutableSupport::EXIT_BAD_ARGUMENTS;
        return exit_code;
    }

    /***********************
     * SIMULATOR PARAMETERS
     ***********************/
    std::string directoryString, filenameString;
    int outputMode; //0 = counts; 1 = mitotic mode events; 2 = mitotic mode sequence sampling
    bool debugOutput;
    unsigned startSeed, endSeed;
    double endTime;
    double normalMu, normalSigma, pPP, pPD, pBC, pAC, pMG; //stochastic model parameters

    //PARSE ARGUMENTS
    directoryString = argv[1];
    filenameString = argv[2];
    outputMode = std::stoi(argv[3]);
    debugOutput = std::stoul(argv[4]);
    startSeed = std::stoul(argv[5]);
    endSeed = std::stoul(argv[6]);
    endTime = std::stod(argv[7]);
    normalMu = std::stod(argv[8]);
    normalSigma = std::stod(argv[9]);
    pPP = std::stod(argv[10]);
    pPD = std::stod(argv[11]);
    pBC = std::stod(argv[12]);
    pAC = std::stod(argv[13]);
    pMG = std::stod(argv[14]);

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

    if (endSeed < startSeed)
    {
        ExecutableSupport::PrintError("Bad start & end seeds (arguments, 5, 6). endSeed must not be < startSeed");
        sane = 0;
    }

    if (endTime <= 0)
    {
        ExecutableSupport::PrintError("Bad endTime (argument 7). endTime must be > 0");
        sane = 0;
    }

    if (normalMu <= 0 || normalSigma <= 0)
    {
        ExecutableSupport::PrintError("Bad cell cycle normal mean or std (arguments 8, 9). Must be  >0");
        sane = 0;
    }

    if (pPP + pPD > 1 || pPP > 1 || pPP < 0 || pPD > 1 || pPD < 0)
    {
        ExecutableSupport::PrintError(
                "Bad mitotic mode probabilities (arguments 10, 11). pPP + pPD should be >=0, <=1, sum should not exceed 1");
        sane = 0;
    }

    if (pBC + pAC + pMG > 1 || pBC > 1 || pBC < 0 || pAC > 1 || pAC < 0 || pMG > 1 || pMG < 0)
    {
        ExecutableSupport::PrintError(
                "Bad specification probabilities (arguments 12, 13, 14). pBC, pAC, pMG should be >=0, <=1, sum should not exceed 1");
        sane = 0;
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
    if (outputMode == 0) *p_log << "Entry\tSeed\tCount\n";
    if (outputMode == 1) *p_log << "Time (hpf)\tSeed\tCellID\tMitotic Mode (0=PP;1=PD;2=DD)\n";
    if (outputMode == 2) *p_log << "Entry\tSeed\tSequence\n";

//Instance RNG
    RandomNumberGenerator* p_RNG = RandomNumberGenerator::Instance();

//Initialise pointers to relevant singleton ProliferativeTypes and Properties
    MAKE_PTR(WildTypeCellMutationState, p_state);
    MAKE_PTR(TransitCellProliferativeType, p_Mitotic);
    MAKE_PTR(DifferentiatedCellProliferativeType, p_PostMitotic);
    MAKE_PTR(RodPhotoreceptor, p_RPh_fate);
    MAKE_PTR(AmacrineCell, p_AC_fate);
    MAKE_PTR(BipolarCell, p_BC_fate);
    MAKE_PTR(MullerGlia, p_MG_fate);
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
        GomesCellCycleModel* p_cycle_model = new GomesCellCycleModel;

        if (debugOutput)
        {
            //Pass ColumnDataWriter to cell cycle model for debug output
            boost::shared_ptr<ColumnDataWriter> p_debugWriter(
                    new ColumnDataWriter(directoryString, filenameString + "DEBUG_" + std::to_string(seed), false, 10));
            p_cycle_model->EnableModelDebugOutput(p_debugWriter);
            debugWriter = &*p_debugWriter;
        }

        //Setup lineages' cycle model with appropriate parameters
        p_cycle_model->SetDimension(2);
        p_cycle_model->SetPostMitoticType(p_PostMitotic);

        //Setup vector containing lineage founder with the properly set up cell cycle model
        std::vector<CellPtr> cells;
        CellPtr p_cell(new Cell(p_state, p_cycle_model));
        p_cell->SetCellProliferativeType(p_Mitotic);
        p_cycle_model->SetModelParameters(normalMu, normalSigma, pPP, pPD, pBC, pAC, pMG);
        p_cycle_model->SetModelProperties(p_RPh_fate, p_AC_fate, p_BC_fate, p_MG_fate);
        if (outputMode == 2) p_cycle_model->EnableSequenceSampler(p_label);
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
        p_simulator->SetDt(0.25);
        p_simulator->SetEndTime(endTime);
        p_simulator->SetOutputDirectory("UnusedSimOutput" + filenameString); //unused output
        p_simulator->Solve();

        //Count lineage size
        unsigned count = cell_population->GetNumRealCells();

        if (outputMode == 0) *p_log << entry_number << "\t" << seed << "\t" << count << "\n";
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
