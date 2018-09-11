#include "../../../../projects/ISP/src/GomesRetinalNeuralFates.hpp"

RodPhotoreceptor::RodPhotoreceptor(unsigned colour)
    : AbstractCellProperty(),
      mColour(colour)
{
}

RodPhotoreceptor::~RodPhotoreceptor()
{
}

unsigned RodPhotoreceptor::GetColour() const
{
    return mColour;
}

AmacrineCell::AmacrineCell(unsigned colour)
    : AbstractCellProperty(),
      mColour(colour)
{
}

AmacrineCell::~AmacrineCell()
{
}

unsigned AmacrineCell::GetColour() const
{
    return mColour;
}

BipolarCell::BipolarCell(unsigned colour)
    : AbstractCellProperty(),
      mColour(colour)
{
}

BipolarCell::~BipolarCell()
{
}

unsigned BipolarCell::GetColour() const
{
    return mColour;
}

MullerGlia::MullerGlia(unsigned colour)
    : AbstractCellProperty(),
      mColour(colour)
{
}

MullerGlia::~MullerGlia()
{
}

unsigned MullerGlia::GetColour() const
{
    return mColour;
}

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(RodPhotoreceptor)
CHASTE_CLASS_EXPORT(AmacrineCell)
CHASTE_CLASS_EXPORT(BipolarCell)
CHASTE_CLASS_EXPORT(MullerGlia)
