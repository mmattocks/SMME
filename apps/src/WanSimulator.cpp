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

#include "CellProliferativeTypesCountWriter.hpp"

int main(int argc, char *argv[])
{
    ExecutableSupport::StartupWithoutShowingCopyright(&argc, &argv);
    //main() returns code indicating sim run success or failure mode
    int exit_code = ExecutableSupport::EXIT_OK;

    if (argc != 23)
    {
        ExecutableSupport::PrintError(
                "Wrong arguments for simulator.\nUsage (replace<> with values, pass bools as 0 or 1):\n WanSimulator <directoryString> <startSeedUnsigned> <endSeedUnsigned> <cmzResidencyTimeDoubleHours> <stemDivisorDouble> <meanProgenitorPopualtion@3dpfDouble> <stdProgenitorPopulation@3dpfDouble> <stemGammaShiftDouble> <stemGammaShapeDouble> <stemGammaScaleDouble> <progenitorGammaShiftDouble> <progenitorGammaShapeDouble> <progenitorGammaScaleDouble> <progenitorSisterShiftDouble> <mMitoticModePhase2Double> <mMitoticModePhase3Double> <pPP1Double(0-1)> <pPD1Double(0-1)> <pPP1Double(0-1)> <pPD1Double(0-1)> <pPP1Double(0-1)> <pPD1Double(0-1)>",
                true);
        exit_code = ExecutableSupport::EXIT_BAD_ARGUMENTS;
        return exit_code;
    }

    /***********************
     * SIMULATOR PARAMETERS
     ***********************/
    std::string directoryString;
    unsigned startSeed, endSeed;
    double cmzResidencyTime;
    double stemDivisor, progenitorMean, progenitorStd;
    double stemGammaShift, stemGammaShape, stemGammaScale, progenitorGammaShift, progenitorGammaShape,
            progenitorGammaScale, progenitorGammaSister;
    double mitoticModePhase2, mitoticModePhase3, pPP1, pPD1, pPP2, pPD2, pPP3, pPD3; //stochastic He model parameters

    //PARSE ARGUMENTS
    directoryString = argv[1];
    startSeed = std::stoul(argv[2]);
    endSeed = std::stoul(argv[3]);
    cmzResidencyTime = std::stod(argv[4]);
    //starting cell number distributions
    stemDivisor = std::stod(argv[5]);
    progenitorMean = std::stod(argv[6]);
    progenitorStd = std::stod(argv[7]);
    //cycle duration params
    stemGammaShift = std::stod(argv[8]);
    stemGammaShape = std::stod(argv[9]);
    stemGammaScale = std::stod(argv[10]);
    progenitorGammaShift = std::stod(argv[11]);
    progenitorGammaShape = std::stod(argv[12]);
    progenitorGammaScale = std::stod(argv[13]);
    progenitorGammaSister = std::stod(argv[14]);
    //He model params
    mitoticModePhase2 = std::stod(argv[15]);
    mitoticModePhase3 = std::stod(argv[16]);
    pPP1 = std::stod(argv[17]);
    pPD1 = std::stod(argv[18]);
    pPP2 = std::stod(argv[19]);
    pPD2 = std::stod(argv[20]);
    pPP3 = std::stod(argv[21]);
    pPD3 = std::stod(argv[22]);

    std::vector<double> stemOffspringParams = { mitoticModePhase2, mitoticModePhase2 + mitoticModePhase3, pPP1, pPD1,
                                                pPP2, pPD2, pPP3, pPD3, progenitorGammaShift, progenitorGammaShape,
                                                progenitorGammaScale, progenitorGammaSister };

    /************************
     * PARAMETER/ARGUMENT SANITY CHECK
     ************************/
    bool sane = 1;

    if (endSeed < startSeed)
    {
        ExecutableSupport::PrintError("Bad start & end seeds (arguments, 3, 4). endSeed must not be < startSeed");
        sane = 0;
    }

    if (cmzResidencyTime <= 0)
    {
        ExecutableSupport::PrintError("Bad CMZ residency time (argument 5). cmzResidencyTime must be positive-valued");
        sane = 0;
    }

    if (stemDivisor <= 0)
    {
        ExecutableSupport::PrintError("Bad stemDivisor (argument 6). stemDivisor must be positive-valued");
        sane = 0;
    }

    if (progenitorMean <= 0 || progenitorStd <= 0)
    {
        ExecutableSupport::PrintError("Bad progenitorMean or progenitorStd (arguments 7,8). Must be positive-valued");
        sane = 0;
    }

    if (stemGammaShift < 0 || stemGammaShape <= 0 || stemGammaScale <= 0)
    {
        ExecutableSupport::PrintError(
                "Bad stemGammaShift, stemGammaShape, or stemGammaScale (arguments 9, 10, 11). Shifts must be >=0, cycle shape and scale params must be positive-valued");
        sane = 0;
    }

    if (progenitorGammaShift < 0 || progenitorGammaShape <= 0 || progenitorGammaScale <= 0 || progenitorGammaSister < 0)
    {
        ExecutableSupport::PrintError(
                "Bad progenitorGammaShift, progenitorGammaShape, progenitorGammaScale, or progenitorGammaSisterShift (arguments 12,13,14,15). Shifts must be >=0, cycle shape and scale params must be positive-valued");
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
     * SIMULATOR SETUP & RUN
     ************************/

    ExecutableSupport::Print("Simulator writing files to directory " + directoryString);

//Instance RNG
    RandomNumberGenerator* p_RNG = RandomNumberGenerator::Instance();

//Initialise pointers to relevant singleton ProliferativeTypes and Properties
    boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
    boost::shared_ptr<AbstractCellProperty> p_Stem(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());
    boost::shared_ptr<AbstractCellProperty> p_Transit(
            CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>());
    boost::shared_ptr<AbstractCellProperty> p_PostMitotic(
            CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>());

//iterate through supplied seed range, executing one simulation per seed
    for (unsigned seed = startSeed; seed <= endSeed; seed++)
    {
        //initialise SimulationTime (permits cellcyclemodel setup)
        SimulationTime::Instance()->SetStartTime(0.0);

        //Reseed the RNG with the required seed
        p_RNG->Reseed(seed);

        //unsigned numberStem = int(std::round(p_RNG->NormalRandomDeviate(stemMean, stemStd)));
        unsigned numberProgenitors = int(std::round(p_RNG->NormalRandomDeviate(progenitorMean, progenitorStd)));
        unsigned numberStem = int(std::round(numberProgenitors / stemDivisor));

        std::vector<CellPtr> stems;
        std::vector<CellPtr> cells;

        for (unsigned i = 0; i < numberStem; i++)
        {
            WanStemCellCycleModel* p_stem_model = new WanStemCellCycleModel;
            p_stem_model->SetDimension(2);
            p_stem_model->SetModelParameters(stemGammaShift, stemGammaShape, stemGammaScale, stemOffspringParams);

            CellPtr p_cell(new Cell(p_state, p_stem_model));
            p_cell->InitialiseCellCycleModel();
            stems.push_back(p_cell);
            cells.push_back(p_cell);
        }

        for (unsigned i = 0; i < numberProgenitors; i++)
        {
            double currTiL = p_RNG->ranf() * cmzResidencyTime;

            HeCellCycleModel* p_prog_model = new HeCellCycleModel;
            p_prog_model->SetDimension(2);
            p_prog_model->SetModelParameters(currTiL, mitoticModePhase2, mitoticModePhase2 + mitoticModePhase3, pPP1,
                                             pPD1, pPP2, pPD2, pPP3, pPD3);
            p_prog_model->EnableKillSpecified();

            CellPtr p_cell(new Cell(p_state, p_prog_model));
            p_cell->InitialiseCellCycleModel();
            cells.push_back(p_cell);
        }

        //Generate 1x#cells mesh for abstract colony
        HoneycombMeshGenerator generator(1, (numberProgenitors + numberStem));
        MutableMesh<2, 2>* p_generating_mesh = generator.GetMesh();
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);

        //Setup cell population
        boost::shared_ptr<NodeBasedCellPopulation<2>> cell_population(new NodeBasedCellPopulation<2>(mesh, cells));
        cell_population->AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();

        //Give Wan stem cells the population & base stem pop size
        for (auto p_cell : stems)
        {
            WanStemCellCycleModel* p_cycle_model = dynamic_cast<WanStemCellCycleModel*>(p_cell->GetCellCycleModel());
            p_cycle_model->EnableExpandingStemPopulation(numberStem, cell_population);
        }

        //Setup simulator & run simulation
        boost::shared_ptr<OffLatticeSimulationPropertyStop<2>> p_simulator(
                new OffLatticeSimulationPropertyStop<2>(*cell_population));
        p_simulator->SetStopProperty(p_Transit); //simulation to stop if no RPCs are left
        p_simulator->SetDt(1);
        p_simulator->SetOutputDirectory(directoryString + "/Seed" + std::to_string(seed) + "Results");
        p_simulator->SetEndTime(8568); // 360dpf - 3dpf simulation start time
        p_simulator->Solve();

        //Reset for next simulation
        SimulationTime::Destroy();
        cell_population.reset();
    }

    p_RNG->Destroy();

    return exit_code;
}
;
