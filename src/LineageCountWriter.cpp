#include "../../../../projects/ISP/src/LineageCountWriter.hpp"

#include "AbstractCellPopulation.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
LineageCountWriter<ELEMENT_DIM, SPACE_DIM>::LineageCountWriter()
    : AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM>("noseed.csv"),
      wSim(),
      wInductionTime(0),
      wSeed(404) // defaults to 404 if no seed provided
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void LineageCountWriter<ELEMENT_DIM, SPACE_DIM>::SetLineageParameters(
        boost::shared_ptr<OffLatticeSimulationPropertyStop<ELEMENT_DIM>> sim, double inductionTime,
        unsigned seed)
{
    wSim = sim;
    wInductionTime = inductionTime;
    wSeed = seed;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void LineageCountWriter<ELEMENT_DIM, SPACE_DIM>::WriteHeader(
        AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    *this->mpOutStream << "I.Time,Seed,Count";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void LineageCountWriter<ELEMENT_DIM, SPACE_DIM>::VisitAnyPopulation(
        AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    SimulationTime* p_simulation_time = SimulationTime::Instance();

    //if the simulation has completed, write the output line
    if (p_simulation_time->IsFinished() || wSim->HasStoppingEventOccurred())
    {
        unsigned count = pCellPopulation->GetNumRealCells();
        *this->mpOutStream << wInductionTime << "," << wSeed << "," << count;
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void LineageCountWriter<ELEMENT_DIM, SPACE_DIM>::Visit(MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{

}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void LineageCountWriter<ELEMENT_DIM, SPACE_DIM>::Visit(CaBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{

}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void LineageCountWriter<ELEMENT_DIM, SPACE_DIM>::Visit(NodeBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{

}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void LineageCountWriter<ELEMENT_DIM, SPACE_DIM>::Visit(PottsBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{

}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void LineageCountWriter<ELEMENT_DIM, SPACE_DIM>::Visit(VertexBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{

}

// Explicit instantiation
template class LineageCountWriter<1, 1> ;
template class LineageCountWriter<1, 2> ;
template class LineageCountWriter<2, 2> ;
template class LineageCountWriter<1, 3> ;
template class LineageCountWriter<2, 3> ;
template class LineageCountWriter<3, 3> ;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(LineageCountWriter)
