#include <cxxtest/TestSuite.h>
#include "../../projects/ISP/src/BoijeCellCycleModel.hpp"
#include "../../projects/ISP/src/BoijeRetinalNeuralFates.hpp"
#include "../../projects/ISP/src/CellPropertiesWriter.hpp"
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
#include "NodeBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"

#include "GeneralisedLinearSpringForce.hpp"
#include "SphereGeometryBoundaryCondition.hpp"
#include "SmartPointers.hpp"



class TestBoijeModel : public AbstractCellBasedTestSuite
{
public:
    void TestBoijeModelFixture()
    {
        EXIT_IF_PARALLEL;

        /*Generate 1x1 mesh for single-cell colony*/
        HoneycombMeshGenerator generator(1, 1);
        MutableMesh < 2, 2 > *p_generating_mesh = generator.GetMesh();

        /*Prepare nodes-only mesh*/
        NodesOnlyMesh < 2 > mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);

        /*Initialise pointers to relevant ProliferativeTypes and Properties*/
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_Mitotic);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_PostMitotic);
        MAKE_PTR(RetinalGanglion, p_RGC_fate);
        MAKE_PTR(AmacrineHorizontal, p_AC_HC_fate);
        MAKE_PTR(ReceptorBipolar, p_PR_BC_fate);

        std::vector < boost::shared_ptr < AbstractCellProperty >> property_write_list;
        property_write_list.push_back(p_Mitotic);
        property_write_list.push_back(p_PostMitotic);
        property_write_list.push_back(p_RGC_fate);
        property_write_list.push_back(p_AC_HC_fate);
        property_write_list.push_back(p_PR_BC_fate);

        BoijeCellCycleModel* p_cycle_model = new BoijeCellCycleModel();
        p_cycle_model->SetDimension(2);
        p_cycle_model->SetPostMitoticType(p_PostMitotic);
        p_cycle_model->SetRGCType(p_RGC_fate);
        p_cycle_model->SetACHCType(p_AC_HC_fate);
        p_cycle_model->SetPRBCType(p_PR_BC_fate);
        p_cycle_model->SetModelParameters(0.32, 0.30, 0.80);

        /*Setup vector containing lineage founder*/
        std::vector < CellPtr > cells;

        CellPtr p_cell(new Cell(p_state, p_cycle_model));
        p_cell->SetCellProliferativeType(p_Mitotic);
        p_cell->SetBirthTime(0.0);
        p_cell->InitialiseCellCycleModel();
        cells.push_back(p_cell);

        NodeBasedCellPopulation < 2 > cell_population(mesh, cells);

        boost::shared_ptr<CellPropertiesWriter<2, 2>> sp_writer(new CellPropertiesWriter<2, 2>());
        sp_writer->SetPropertiesToWrite(property_write_list);
        cell_population.AddPopulationWriter(sp_writer);

        OffLatticeSimulationPropertyStop < 2 > simulator(cell_population);
        simulator.SetStopProperty(p_Mitotic);
        simulator.SetOutputDirectory("BoijeProto");
        simulator.SetSamplingTimestepMultiple(120);
        simulator.SetEndTime(30.0);

        /*Pass force law to simulator*/
        //MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        //simulator.AddForce(p_force);
        simulator.Solve();

    }

};
