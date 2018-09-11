#ifndef BOIJERETINALNEURALFATES_HPP_
#define BOIJERETINALNEURALFATES_HPP_

#include <boost/shared_ptr.hpp>
#include "AbstractCellProperty.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

class RetinalGanglion : public AbstractCellProperty
{
private:

    /**
     * Colour for use by visualizer.
     */
    unsigned mColour;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellProperty>(*this);
        archive & mColour;
    }

public:

    /**
     * Constructor.
     *
     * @param colour  what colour cells with this property should be in the visualizer (defaults to 6)
     */
    RetinalGanglion(unsigned colour=3);

    /**
     * Destructor.
     */
    virtual ~RetinalGanglion();

    /**
     * @return #mColour.
     */
    unsigned GetColour() const;
};

class AmacrineHorizontal : public AbstractCellProperty
{
private:

    /**
     * Colour for use by visualizer.
     */
    unsigned mColour;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellProperty>(*this);
        archive & mColour;
    }

public:

    /**
     * Constructor.
     *
     * @param colour  what colour cells with this property should be in the visualizer (defaults to 6)
     */
    AmacrineHorizontal(unsigned colour=4);

    /**
     * Destructor.
     */
    virtual ~AmacrineHorizontal();

    /**
     * @return #mColour.
     */
    unsigned GetColour() const;
};

class ReceptorBipolar : public AbstractCellProperty
{
private:

    /**
     * Colour for use by visualizer.
     */
    unsigned mColour;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellProperty>(*this);
        archive & mColour;
    }

public:

    /**
     * Constructor.
     *
     * @param colour  what colour cells with this property should be in the visualizer (defaults to 6)
     */
    ReceptorBipolar(unsigned colour=5);

    /**
     * Destructor.
     */
    virtual ~ReceptorBipolar();

    /**
     * @return #mColour.
     */
    unsigned GetColour() const;
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(RetinalGanglion)
CHASTE_CLASS_EXPORT(AmacrineHorizontal)
CHASTE_CLASS_EXPORT(ReceptorBipolar)

#endif /* BOIJERETINALNEURALFATES_HPP_ */
