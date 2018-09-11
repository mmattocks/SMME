#include <cxxtest/TestSuite.h>
#include "../../projects/ISP/src/GomesCellCycleModel.hpp"
#include "../../projects/ISP/src/GomesRetinalNeuralFates.hpp"
#include "../../projects/ISP/src/OffLatticeSimulationPropertyStop.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "PetscSetupAndFinalize.hpp"

#include "AbstractCellBasedTestSuite.hpp"

#include "HoneycombMeshGenerator.hpp"
#include "NodesOnlyMesh.hpp"
#include "HoneycombVertexMeshGenerator.hpp"

#include "CellsGenerator.hpp"

#include "WildTypeCellMutationState.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "CellLabel.hpp"

#include "NodeBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"

#include "GeneralisedLinearSpringForce.hpp"
#include "SphereGeometryBoundaryCondition.hpp"
#include "SmartPointers.hpp"

#include "ColumnDataWriter.hpp"
#include "LogFile.hpp"

class TestKolmogorov : public AbstractCellBasedTestSuite
{
public:
    bool debugOutput = true;

    /************************
     * SIMULATION PARAMETERS
     ************************/

    //Define start and end RNG seeds; determines:
    //#lineages
    //unique sequence of RNG results for each lineage
    const unsigned startSeed = 283;
    const unsigned endSeed = 284;

    /**********************
     * GOMES SSM PARAMETERS
     **********************/
    //cell cycle duration - parameters are for normal distribution (x), setting cell cycle time exp(x).
    //default params give lognormal PDF with mean 56 hr std 18.6hr
    const double normalMu = 3.9716;
    const double normalSigma = .32839;
    const double PP = .055;
    const double PD = .221;

    void TestGomesKolmogorov()
    {
        Timer::Reset();
        Timer::Print("Begin TestGomesKolmogorov @ ");

        /*Initialise pointers to relevant ProliferativeTypes and Properties*/
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_Mitotic);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_PostMitotic);
        MAKE_PTR(RodPhotoreceptor, p_RPh_fate);
        MAKE_PTR(AmacrineCell, p_AC_fate);
        MAKE_PTR(BipolarCell, p_BC_fate);
        MAKE_PTR(MullerGlia, p_MG_fate);
        MAKE_PTR(CellLabel, p_label);

        //Setup sequence logfile
        LogFile* p_log = LogFile::Instance();
        p_log->Set(0, "TestKolmogorovStrings", "GomesStrings");
        *p_log << "Entry\tSeed\tSequence\n";

        //entry counter for sequence writer
        unsigned entry_number = 1;

        //Setup RNG
        RandomNumberGenerator* p_RNG = RandomNumberGenerator::Instance();

        //iterate through supplied seed range, executing one simulation per seed
        for (unsigned seed = startSeed; seed <= endSeed; seed++)
        {
            //write seed to log - sequence written by cellcyclemodel objects
            *p_log << entry_number << "\t" << seed << "\t";

            //initialise pointer to debugWriter
            ColumnDataWriter* debugWriter;

            //Reseed the RNG with the required seed
            p_RNG->Reseed(seed);

            //Setup a Gomes cell cycle model
            GomesCellCycleModel* p_cycle_model = new GomesCellCycleModel;
            p_cycle_model->SetPostMitoticType(p_PostMitotic);
            p_cycle_model->SetModelProperties(p_RPh_fate, p_AC_fate, p_BC_fate, p_MG_fate);
            p_cycle_model->SetModelParameters(normalMu, normalSigma, PP, PD);
            p_cycle_model->EnableSequenceSampler(p_label);

            if (debugOutput)
            {
                //Setup the writer and model for debug output
                debugWriter = SetupDebugOutput(p_cycle_model, "Gomes", seed);
            }

            /*Setup vector containing lineage founder*/
            std::vector<CellPtr> cells;
            CellPtr p_cell(new Cell(p_state, p_cycle_model));
            p_cell->SetCellProliferativeType(p_Mitotic);
            p_cell->AddCellProperty(p_label);
            p_cell->InitialiseCellCycleModel();
            cells.push_back(p_cell);

            /*Generate 1x1 mesh for single-cell colony*/
            HoneycombMeshGenerator generator(1, 1);
            MutableMesh<2, 2> *p_generating_mesh = generator.GetMesh();
            /*Prepare nodes-only mesh*/
            NodesOnlyMesh<2> mesh;
            mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);

            //Setup cell population
            NodeBasedCellPopulation<2>* cell_population(new NodeBasedCellPopulation<2>(mesh, cells));

            //Setup simulator & run simulation
            boost::shared_ptr<OffLatticeSimulationPropertyStop<2>> p_simulator(
                    new OffLatticeSimulationPropertyStop<2>(*cell_population));
            p_simulator->SetStopProperty(p_Mitotic); //simulation to stop if no mitotic cells are left
            p_simulator->SetDt(0.25);
            p_simulator->SetEndTime(240);
            p_simulator->SetOutputDirectory("UnusedSimOutputTestKolmogorov"); //unused output
            p_simulator->Solve();

            //Write newline to log
            *p_log << "\n";

            //Count lineage size
            unsigned count = cell_population->GetNumRealCells();
            //Output simulation report to console
            std::stringstream reportStream;
            reportStream << "Simulation Seed " << seed << ": Final Lineage Count: " << count
                    << " Sim Complete @ ";
            auto simReport = reportStream.str();
            Timer::Print(simReport);

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

    ColumnDataWriter* SetupDebugOutput(GomesCellCycleModel* p_model, std::string model, unsigned seed)
    {
        //Pass ColumnDataWriter to cell cycle model for debug output
        std::string debugDir = "TestKolmogorovDebugFiles";
        std::stringstream filenameStream;
        filenameStream << "TestKolmogorovDebug" << model << "_S_" << seed;
        std::string debugFilename = filenameStream.str();
        boost::shared_ptr<ColumnDataWriter> p_debugWriter(
                new ColumnDataWriter(debugDir, debugFilename, false, 10));
        p_model->EnableModelDebugOutput(p_debugWriter);
        ColumnDataWriter* debugWriter = &*p_debugWriter;
        return debugWriter;
    }
};
