#ifndef LINEAGECOUNTWRITER_HPP_
#define LINEAGECOUNTWRITER_HPP_

#include "AbstractCellPopulationWriter.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "OffLatticeSimulationPropertyStop.hpp"
#include "SmartPointers.hpp"

/**
 * A class written using the visitor pattern for writing the final count of a cell population,
 * intended for use with He and Boije cell cycle models and the OffLatticeSimulationPropertyStop
 * simulation class
 *
 * The output file is called results.vizcelltypes by default. If VTK is switched on,
 * then the writer also specifies the VTK output for each cell, which is stored in
 * the VTK cell data "Cell types" by default.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class LineageCountWriter : public AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM>
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
        archive & wInductionTime;
        archive & wSeed;
    }
    
protected:

public:
    boost::shared_ptr<OffLatticeSimulationPropertyStop<ELEMENT_DIM>> wSim;
    double wInductionTime;
    unsigned wSeed;

    /**
     * Default constructor.
     */
    LineageCountWriter();

    /**
     * Provide current simulation's induction time and seed variables for output
     */
    void SetLineageParameters(boost::shared_ptr<OffLatticeSimulationPropertyStop<ELEMENT_DIM>> sim, double inductionTime, unsigned seed);

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
EXPORT_TEMPLATE_CLASS_ALL_DIMS(LineageCountWriter)

#endif /* LINEAGECOUNTWRITER_HPP_ */
