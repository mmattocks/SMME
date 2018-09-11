#include "../../../../projects/ISP/src/HeAth5Mo.hpp"

Ath5Mo::Ath5Mo(unsigned colour)
    : AbstractCellProperty(),
      mColour(colour)
{
}

Ath5Mo::~Ath5Mo()
{
}

unsigned Ath5Mo::GetColour() const
{
    return mColour;
}

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(Ath5Mo)
