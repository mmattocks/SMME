#include <iostream>
#include <string>

#include <cxxtest/TestSuite.h>
#include "ExecutableSupport.hpp"
#include "Exception.hpp"
#include "PetscTools.hpp"
#include "PetscException.hpp"

#include "WanStemCellCycleModel.hpp"
#include "HeCellCycleModel.hpp"
#include "OffLatticeSimulationPropertyStop.hpp"

#include "AbstractCellBasedTestSuite.hpp"

#include "WildTypeCellMutationState.hpp"
#include "StemCellProliferativeType.hpp"
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

    if (argc != 32)
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
    bool debugOutput;
    unsigned startSeed, endSeed;
    double endTime, cmzResidencyTime;
    double stemMean, stemStd, progenitorMean, progenitorStd;
    double stemRatePeak, stemIncreasingSlope, stemDecreasingSlope, progenitorRatePeak, progenitorIncreasingSlope,
            progenitorDecreasingSlope;
    double stemGammaShift, stemGammaShape, stemGammaScale, progenitorGammaShift, progenitorGammaShape,
            progenitorGammaScale, progenitorGammaSister;
    double mitoticModePhase2, mitoticModePhase3, pPP1, pPD1, pPP2, pPD2, pPP3, pPD3; //stochastic He model parameters

    //PARSE ARGUMENTS
    directoryString = argv[1];
    filenameString = argv[2];
    debugOutput = std::stoul(argv[3]);
    startSeed = std::stoul(argv[4]);
    endSeed = std::stoul(argv[5]);
    endTime = std::stod(argv[6]);
    cmzResidencyTime = std::stod(argv[7]);
    //starting cell number distributions
    stemMean = std::stod(argv[8]);
    stemStd = std::stod(argv[9]);
    progenitorMean = std::stod(argv[10]);
    progenitorStd = std::stod(argv[11]);
    //adjustable cycle duration params
    stemRatePeak = std::stod(argv[12]);
    stemIncreasingSlope = std::stod(argv[13]);
    stemDecreasingSlope = std::stod(argv[14]);
    progenitorRatePeak = std::stod(argv[15]);
    progenitorIncreasingSlope = std::stod(argv[16]);
    progenitorDecreasingSlope = std::stod(argv[17]);
    //cycle duration params
    stemGammaShift = std::stod(argv[18]);
    stemGammaShape = std::stod(argv[19]);
    stemGammaScale = std::stod(argv[20]);
    progenitorGammaShift = std::stod(argv[21]);
    progenitorGammaShape = std::stod(argv[22]);
    progenitorGammaScale = std::stod(argv[23]);
    progenitorGammaSister = std::stod(argv[24]);
    //He model params
    mitoticModePhase2 = std::stod(argv[25]);
    mitoticModePhase3 = std::stod(argv[26]);
    pPP1 = std::stod(argv[27]);
    pPD1 = std::stod(argv[28]);
    pPP2 = std::stod(argv[29]);
    pPD2 = std::stod(argv[30]);
    pPP3 = std::stod(argv[31]);
    pPD3 = std::stod(argv[32]);

    std::vector<double> stemOffspringParams = { mitoticModePhase2, mitoticModePhase2+mitoticModePhase3, pPP1, pPD1, pPP2, pPD2, pPP3, pPD3, progenitorGammaShift, progenitorGammaShape, progenitorGammaScale, progenitorGammaSister};


    /************************
     * PARAMETER/ARGUMENT SANITY CHECK
     ************************/
    bool sane = 1;

    if (endSeed < startSeed)
    {
        ExecutableSupport::PrintError("Bad start & end seeds (arguments, 8, 9). endSeed must not be < startSeed");
        sane = 0;
    }

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

//Write appropriate header to log
*p_log << "Entry\tSeed\tTotalCells\tStemCount\tProgenitorCount\tPostMitoticCount\n";

//Instance RNG
RandomNumberGenerator* p_RNG = RandomNumberGenerator::Instance();

//Initialise pointers to relevant singleton ProliferativeTypes and Properties
MAKE_PTR(WildTypeCellMutationState, p_state);
MAKE_PTR(StemCellProliferativeType, p_Stem);
MAKE_PTR(TransitCellProliferativeType, p_Transit);
MAKE_PTR(DifferentiatedCellProliferativeType, p_PostMitotic);
MAKE_PTR(Ath5Mo, p_Morpholino);
MAKE_PTR(CellLabel, p_label);

/************************
 * SIMULATOR SETUP & RUN
 ************************/

//iterate through supplied seed range, executing one simulation per seed
for (unsigned seed = startSeed; seed <= endSeed; seed++)
{
    //initialise pointer to debugWriter
    boost::shared_ptr<ColumnDataWriter> p_debugWriter(
            new ColumnDataWriter(directoryString, filenameString + "DEBUG_" + std::to_string(seed), false, 10));

    //initialise SimulationTime (permits cellcyclemodel setup)
    SimulationTime::Instance()->SetStartTime(0.0);

    //Reseed the RNG with the required seed
    p_RNG->Reseed(seed);

    unsigned numberStem = int(std::round(p_RNG->NormalRandomDeviate(stemMean, stemStd)));
    unsigned numberProgenitors = int(std::round(p_RNG->NormalRandomDeviate(progenitorMean, progenitorStd)));

    std::vector<CellPtr> cells;

    for(unsigned i=0;i<numberStem;i++)
    {
        WanStemCellCycleModel* p_stem_model = new WanStemCellCycleModel;

        if (debugOutput)
        {
            p_stem_model->EnableModelDebugOutput(p_debugWriter);
        }

        p_stem_model->SetDimension(2);
        p_stem_model->SetTransitType(p_Transit);
        p_stem_model->SetModelParameters(stemGammaShift, stemGammaShape, stemGammaScale, stemOffspringParams);
        p_stem_model->SetTimeDependentCycleDuration(stemRatePeak, stemIncreasingSlope, stemDecreasingSlope);

        CellPtr p_cell(new Cell(p_state, p_stem_model));
        p_cell->SetCellProliferativeType(p_Stem);
        p_cell->InitialiseCellCycleModel();
        cells.push_back(p_cell);
    }

    for(unsigned i=0;i<numberProgenitors;i++)
    {
        HeCellCycleModel* p_prog_model = new HeCellCycleModel;

        if (debugOutput)
        {
            p_prog_model->EnableModelDebugOutput(p_debugWriter);
        }

        double currTiL = p_RNG->ranf() * cmzResidencyTime;

        p_prog_model->SetDimension(2);
        p_prog_model->SetPostMitoticType(p_PostMitotic);
        p_prog_model->SetModelParameters(currTiL, mitoticModePhase2, mitoticModePhase2 + mitoticModePhase3, pPP1,
                pPD1, pPP2, pPD2, pPP3, pPD3);
        p_prog_model->SetTimeDependentCycleDuration(progenitorRatePeak, progenitorIncreasingSlope, progenitorDecreasingSlope);

        CellPtr p_cell(new Cell(p_state, p_prog_model));
        p_cell->SetCellProliferativeType(p_Transit);
        p_cell->InitialiseCellCycleModel();
        cells.push_back(p_cell);
    }

    //Generate 1x1 mesh for abstract colony
    HoneycombMeshGenerator generator(1, 1);
    MutableMesh<2, 2>* p_generating_mesh = generator.GetMesh();
    NodesOnlyMesh<2> mesh;
    mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);

    //Setup cell population
    NodeBasedCellPopulation<2>* cell_population(new NodeBasedCellPopulation<2>(mesh, cells));

    //Setup simulator & run simulation
    boost::shared_ptr<OffLatticeSimulationPropertyStop<2>> p_simulator(
            new OffLatticeSimulationPropertyStop<2>(*cell_population));
    p_simulator->SetStopProperty(p_Transit);//simulation to stop if no RPCs are left
    p_simulator->SetDt(0.05);
    p_simulator->SetEndTime(endTime);
    p_simulator->SetOutputDirectory("UnusedSimOutput" + filenameString);//unused output
    p_simulator->Solve();

    //Count population size and composition
    unsigned count = cell_population->GetNumRealCells();
    unsigned stem = p_Stem->GetCellCount();
    unsigned progenitor = p_Transit->GetCellCount();
    unsigned postmitotic = p_PostMitotic->GetCellCount();

    *p_log << entry_number << "\t" << seed << "\t" << count << "\t" << stem << "\t" << progenitor << "\t" << postmitotic << "\n";

    //Reset for next simulation
    SimulationTime::Destroy();
    delete cell_population;
    entry_number++;

    if (debugOutput)
    {
        p_debugWriter->Close();
    }

}

p_RNG->Destroy();
LogFile::Close();

return exit_code;
}
;
