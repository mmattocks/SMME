#ifndef GOMESRETINALNEURALFATES_HPP_
#define GOMESRETINALNEURALFATES_HPP_

#include <boost/shared_ptr.hpp>
#include "AbstractCellProperty.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

class RodPhotoreceptor : public AbstractCellProperty
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
    RodPhotoreceptor(unsigned colour=3);

    /**
     * Destructor.
     */
    virtual ~RodPhotoreceptor();

    /**
     * @return #mColour.
     */
    unsigned GetColour() const;
};

class AmacrineCell : public AbstractCellProperty
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
    AmacrineCell(unsigned colour=4);

    /**
     * Destructor.
     */
    virtual ~AmacrineCell();

    /**
     * @return #mColour.
     */
    unsigned GetColour() const;
};

class BipolarCell : public AbstractCellProperty
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
    BipolarCell(unsigned colour=5);

    /**
     * Destructor.
     */
    virtual ~BipolarCell();

    /**
     * @return #mColour.
     */
    unsigned GetColour() const;
};

class MullerGlia : public AbstractCellProperty
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
    MullerGlia(unsigned colour=6);

    /**
     * Destructor.
     */
    virtual ~MullerGlia();

    /**
     * @return #mColour.
     */
    unsigned GetColour() const;
};
#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(RodPhotoreceptor)
CHASTE_CLASS_EXPORT(AmacrineCell)
CHASTE_CLASS_EXPORT(BipolarCell)
CHASTE_CLASS_EXPORT(MullerGlia)

#endif /* GOMESRETINALNEURALFATES_HPP_ */
