#include "../../../../projects/ISP/src/BoijeRetinalNeuralFates.hpp"

RetinalGanglion::RetinalGanglion(unsigned colour)
    : AbstractCellProperty(),
      mColour(colour)
{
}

RetinalGanglion::~RetinalGanglion()
{
}

unsigned RetinalGanglion::GetColour() const
{
    return mColour;
}

AmacrineHorizontal::AmacrineHorizontal(unsigned colour)
    : AbstractCellProperty(),
      mColour(colour)
{
}

AmacrineHorizontal::~AmacrineHorizontal()
{
}

unsigned AmacrineHorizontal::GetColour() const
{
    return mColour;
}

ReceptorBipolar::ReceptorBipolar(unsigned colour)
    : AbstractCellProperty(),
      mColour(colour)
{
}

ReceptorBipolar::~ReceptorBipolar()
{
}

unsigned ReceptorBipolar::GetColour() const
{
    return mColour;
}

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(RetinalGanglion)
CHASTE_CLASS_EXPORT(AmacrineHorizontal)
CHASTE_CLASS_EXPORT(ReceptorBipolar)
