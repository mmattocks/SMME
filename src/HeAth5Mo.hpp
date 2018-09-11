#ifndef HEATH5MO_HPP_
#define HEATH5MO_HPP_

#include <boost/shared_ptr.hpp>
#include "AbstractCellProperty.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

class Ath5Mo : public AbstractCellProperty
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
    Ath5Mo(unsigned colour=3);

    /**
     * Destructor.
     */
    virtual ~Ath5Mo();

    /**
     * @return #mColour.
     */
    unsigned GetColour() const;
};


#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(Ath5Mo)
#endif /* HEATH5MO_HPP */
