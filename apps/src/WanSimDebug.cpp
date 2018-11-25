#include <iostream>
#include <string>

#include <cxxtest/TestSuite.h>
#include "ExecutableSupport.hpp"
#include "Exception.hpp"
#include "PetscTools.hpp"
#include "PetscException.hpp"

#include "WanStemCellCycleModel.hpp"
#include "HeCellCycleModel.hpp"
#include "OffLatticeSimulation.hpp"

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

    /***********************
     * SIMULATOR PARAMETERS
     ***********************/
    std::string directoryString, filenameString;

    //PARSE ARGUMENTS
    directoryString = argv[1];
    filenameString = argv[2];

    /************************
     * SIMULATOR OUTPUT SETUP
     ************************/

//Set up singleton LogFile
    LogFile* p_log = LogFile::Instance();
    p_log->Set(0, directoryString, filenameString);

    ExecutableSupport::Print("Simulator writing file " + filenameString + " to directory " + directoryString);

//Log entry counter
    //unsigned entry_number = 1;

//Write appropriate header to log
    *p_log << "Entry\tTotalCells\tStemCount\tProgenitorCount\tPostMitoticCount\n";

//Instance RNG
    RandomNumberGenerator* p_RNG = RandomNumberGenerator::Instance();

//Initialise pointers to relevant singleton ProliferativeTypes and Properties
    boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
    boost::shared_ptr<AbstractCellProperty> p_Stem(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());
    boost::shared_ptr<AbstractCellProperty> p_Transit(CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>());
    boost::shared_ptr<AbstractCellProperty> p_PostMitotic(CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>());

    /************************
     * SIMULATOR SETUP & RUN
     ************************/

    //initialise SimulationTime (permits cellcyclemodel setup)
    SimulationTime::Instance()->SetStartTime(0.0);

    std::vector<CellPtr> cells;


    WanStemCellCycleModel* p_stem_model = new WanStemCellCycleModel;

    boost::shared_ptr<ColumnDataWriter> p_debugWriter(
            new ColumnDataWriter(directoryString, filenameString + "DEBUG_WAN", false, 10));

    p_stem_model->SetDimension(2);
    p_stem_model->EnableModelDebugOutput(p_debugWriter);
    CellPtr p_cell(new Cell(p_state, p_stem_model));
    p_cell->InitialiseCellCycleModel();
    cells.push_back(p_cell);

    //Generate 1x#cells mesh for abstract colony
    HoneycombMeshGenerator generator(1, 1);
    MutableMesh<2, 2>* p_generating_mesh = generator.GetMesh();
    NodesOnlyMesh<2> mesh;
    mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);

    //Setup cell population
    boost::shared_ptr<NodeBasedCellPopulation<2>> cell_population(new NodeBasedCellPopulation<2>(mesh, cells));
    cell_population->AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();

    p_stem_model->EnableExpandingStemPopulation(1, cell_population);

    //Setup simulator & run simulation
    boost::shared_ptr<OffLatticeSimulation<2>> p_simulator(
            new OffLatticeSimulation<2>(*cell_population));
    p_simulator->SetDt(1);
    p_simulator->SetOutputDirectory(directoryString + "/WanDebug");
    p_simulator->SetEndTime(100);
    p_simulator->Solve();

    std::vector<unsigned> prolTypes = cell_population->GetCellProliferativeTypeCount();
    unsigned realCells = cell_population->GetNumRealCells();
    unsigned allCells = cell_population->GetNumAllCells();

    Timer::Print("stem: " + std::to_string(prolTypes[0]) + " transit: " + std::to_string(prolTypes[1]) + " postmitotic: " + std::to_string(prolTypes[2]));
    Timer::Print("realcells: " + std::to_string(realCells) + " allcells: " + std::to_string(allCells));

    //Reset for next simulation
    SimulationTime::Destroy();
    cell_population.reset();

    p_RNG->Destroy();
    LogFile::Close();

    return exit_code;
}
;
