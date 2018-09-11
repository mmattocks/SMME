#ifndef CELLPROPERTIESWRITER_HPP_
#define CELLPROPERTIESWRITER_HPP_

#include "AbstractCellPopulationWriter.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "SmartPointers.hpp"
#include <vector>
#include "AbstractCellProperty.hpp"

/**
 * A class written using the visitor pattern for writing cell proliferative types to file.
 *
 * The output file is called results.vizcelltypes by default. If VTK is switched on,
 * then the writer also specifies the VTK output for each cell, which is stored in
 * the VTK cell data "Cell types" by default.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class CellPropertiesWriter : public AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM>
{
private:
/*
    ** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM> >(*this);
        archive & v_properties;
    }
    
protected:

std::vector<boost::shared_ptr<AbstractCellProperty>> v_properties;

public:

    /**
     * Default constructor.
     */
    CellPropertiesWriter();

	void SetPropertiesToWrite(std::vector<boost::shared_ptr<AbstractCellProperty>> properties_to_write);

    /**

     * @param pCellPopulation a pointer to the cell population owning the cell
     */
    void VisitAnyPopulation(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);
    
    /**
     * Write the header to file.
     *
     * @param pCellPopulation a pointer to the population to be written.
     */
    void WriteHeader(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);

    /**
     * <<Overridden/nonfunctional>>Visit the population and write the data.
     * @param pCellPopulation a pointer to the MeshBasedCellPopulation to visit.
     */
    void Visit(MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);

    /**
     * <<Overridden/nonfunctional>>Visit the population and write the data.
     * @param pCellPopulation a pointer to the CaBasedCellPopulation to visit.
     */
    void Visit(CaBasedCellPopulation<SPACE_DIM>* pCellPopulation);

    /**
     *  <<Overridden/nonfunctional>>Visit the population and write the data.
     * @param pCellPopulation a pointer to the NodeBasedCellPopulation to visit.
     */
    void Visit(NodeBasedCellPopulation<SPACE_DIM>* pCellPopulation);

    /**
     *  <<Overridden/nonfunctional>>Visit the population and write the data.
     * @param pCellPopulation a pointer to the PottsBasedCellPopulation to visit.
     */
    void Visit(PottsBasedCellPopulation<SPACE_DIM>* pCellPopulation);

    /**
     *  <<Overridden/nonfunctional>>Visit the population and write the data.
     * @param pCellPopulation a pointer to the VertexBasedCellPopulation to visit.
     */
    void Visit(VertexBasedCellPopulation<SPACE_DIM>* pCellPopulation);
    
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CellPropertiesWriter)

#endif /* CELLPROPERTIESWRITER_HPP_ */
