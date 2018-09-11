#include "../../../../projects/ISP/src/CellPropertiesWriter.hpp"

#include "AbstractCellPopulation.hpp"
#include "AbstractCellProperty.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellPropertiesWriter<ELEMENT_DIM, SPACE_DIM>::CellPropertiesWriter()
    : AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM>("celltypes.csv"),
    v_properties()
    
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellPropertiesWriter<ELEMENT_DIM, SPACE_DIM>::VisitAnyPopulation(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    int i = 0;
    
	for (std::vector<boost::shared_ptr<AbstractCellProperty>>::iterator it = v_properties.begin() ; it != v_properties.end(); it++, i++)
	{
		boost::shared_ptr<AbstractCellProperty> p_property = v_properties[i];
		unsigned count = p_property->GetCellCount();
		*this->mpOutStream << count << ",";
		
	}
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellPropertiesWriter<ELEMENT_DIM, SPACE_DIM>::SetPropertiesToWrite(std::vector<boost::shared_ptr<AbstractCellProperty>> properties_to_write)
{

v_properties = properties_to_write;
	
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellPropertiesWriter<ELEMENT_DIM, SPACE_DIM>::WriteHeader(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
	*this->mpOutStream << "Time,";
	int i = 0;
	for (std::vector<boost::shared_ptr<AbstractCellProperty>>::iterator it = v_properties.begin() ; it != v_properties.end(); it++, i++)
	{
		*this->mpOutStream << typeid(*v_properties[i]).name() << ",";
	}
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellPropertiesWriter<ELEMENT_DIM, SPACE_DIM>::Visit(MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    int i = 0;
    
	for (std::vector<boost::shared_ptr<AbstractCellProperty>>::iterator it = v_properties.begin() ; it != v_properties.end(); it++, i++)
	{
		boost::shared_ptr<AbstractCellProperty> p_property = v_properties[i];
		unsigned count = p_property->GetCellCount();
		*this->mpOutStream << count << ",";
		
	}
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellPropertiesWriter<ELEMENT_DIM, SPACE_DIM>::Visit(CaBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    int i = 0;
    
	for (std::vector<boost::shared_ptr<AbstractCellProperty>>::iterator it = v_properties.begin() ; it != v_properties.end(); it++, i++)
	{
		boost::shared_ptr<AbstractCellProperty> p_property = v_properties[i];
		unsigned count = p_property->GetCellCount();
		*this->mpOutStream << count << ",";
		
	}
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellPropertiesWriter<ELEMENT_DIM, SPACE_DIM>::Visit(NodeBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    int i = 0;
    
	for (std::vector<boost::shared_ptr<AbstractCellProperty>>::iterator it = v_properties.begin() ; it != v_properties.end(); it++, i++)
	{
		boost::shared_ptr<AbstractCellProperty> p_property = v_properties[i];
		unsigned count = p_property->GetCellCount();
		*this->mpOutStream << count << ",";
		
	}
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellPropertiesWriter<ELEMENT_DIM, SPACE_DIM>::Visit(PottsBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    int i = 0;
    
	for (std::vector<boost::shared_ptr<AbstractCellProperty>>::iterator it = v_properties.begin() ; it != v_properties.end(); it++, i++)
	{
		boost::shared_ptr<AbstractCellProperty> p_property = v_properties[i];
		unsigned count = p_property->GetCellCount();
		*this->mpOutStream << count << ",";
		
	}
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellPropertiesWriter<ELEMENT_DIM, SPACE_DIM>::Visit(VertexBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    int i = 0;
    
	for (std::vector<boost::shared_ptr<AbstractCellProperty>>::iterator it = v_properties.begin() ; it != v_properties.end(); it++, i++)
	{
		boost::shared_ptr<AbstractCellProperty> p_property = v_properties[i];
		unsigned count = p_property->GetCellCount();
		*this->mpOutStream << count << ",";
		
	}
}

// Explicit instantiation
template class CellPropertiesWriter<1,1>;
template class CellPropertiesWriter<1,2>;
template class CellPropertiesWriter<2,2>;
template class CellPropertiesWriter<1,3>;
template class CellPropertiesWriter<2,3>;
template class CellPropertiesWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CellPropertiesWriter)
