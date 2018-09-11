#ifndef NEURALTFS_HPP_
#define NEURALTFS_HPP_

#include <boost/shared_ptr.hpp>
#include "AbstractCellProperty.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

class Atoh7 : public AbstractCellProperty
{
private:
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
    }

public:

    /**
     * Constructor.
     */
    Atoh7();

    /**
     * Destructor.
     */
    virtual ~Atoh7();
};

class Ptf1a : public AbstractCellProperty
{
private:

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
    }

public:

    /**
     * Constructor.
     *
     */
    Ptf1a();

    /**
     * Destructor.
     */
    virtual ~Ptf1a();
};

class vsx2 : public AbstractCellProperty
{
private:

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
    }

public:

    /**
     * Constructor.
     *
     */
    vsx2();

    /**
     * Destructor.
     */
    virtual ~vsx2();
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(Atoh7)
CHASTE_CLASS_EXPORT(Ptf1a)
CHASTE_CLASS_EXPORT(vsx2)

#endif /* NEURALTFS_HPP_ */
