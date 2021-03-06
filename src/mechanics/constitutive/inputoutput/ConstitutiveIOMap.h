#pragma once

#include <map>
#include <memory>

#include "mechanics/constitutive/inputoutput/ConstitutiveIOBase.h"

namespace NuTo
{

namespace Constitutive
{
enum class eInput;
enum class eOutput;
}


template <typename IOEnum>
class ConstitutiveIOMap : public std::map<IOEnum, std::unique_ptr<ConstitutiveIOBase>>
{
public:
    ConstitutiveIOMap() = default;
    ConstitutiveIOMap(const ConstitutiveIOMap& other);
    NuTo::ConstitutiveIOMap<IOEnum>& Merge(const ConstitutiveIOMap& other);
    bool Contains(IOEnum rEnum) const;
    template <int TDim>
    void Add(IOEnum rEnum)
    {
        this->operator[](rEnum) = NuTo::ConstitutiveIOBase::makeConstitutiveIO<TDim>(rEnum);
    }
};

using ConstitutiveInputMap = ConstitutiveIOMap<Constitutive::eInput>;
using ConstitutiveOutputMap = ConstitutiveIOMap<Constitutive::eOutput>;
} // namespace NuTo
