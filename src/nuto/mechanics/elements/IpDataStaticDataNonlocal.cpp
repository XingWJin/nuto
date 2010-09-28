// $ld: $ 
// IpDataStaticDataNonlocal.cpp
// created Apr 29, 2010 by Joerg F. Unger

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/ConstitutiveStaticDataBase.h"
#include "nuto/mechanics/elements/IpDataStaticDataNonlocal.h"
#include "nuto/mechanics/elements/ElementBase.h"
NuTo::IpDataStaticDataNonlocal::IpDataStaticDataNonlocal() :NuTo::IpDataBase::IpDataBase() ,
    NuTo::IpDataStaticDataBase::IpDataStaticDataBase() , NuTo::IpDataNonlocalBase::IpDataNonlocalBase()
{

}

NuTo::IpDataStaticDataNonlocal::~IpDataStaticDataNonlocal()
{
	if (mStaticData!=0)
		delete mStaticData;
}

void NuTo::IpDataStaticDataNonlocal::Initialize(const ElementBase* rElement, const ConstitutiveBase* rConstitutive)
{
	if (mStaticData!=0)
		delete mStaticData;
	if (rConstitutive!=0)
	    mStaticData = rElement->AllocateStaticData(rConstitutive);
	else
		mStaticData = 0;
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::IpDataStaticDataNonlocal::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::IpDataStaticDataNonlocal::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::IpDataStaticDataNonlocal::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::IpDataStaticDataNonlocal::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::IpDataStaticDataNonlocal::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::IpDataStaticDataNonlocal::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::IpDataStaticDataNonlocal::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize IpDataStaticDataNonlocal" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(IpDataStaticDataBase)
       & BOOST_SERIALIZATION_BASE_OBJECT_NVP(IpDataNonlocalBase);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize IpDataStaticDataNonlocal" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::IpDataStaticDataNonlocal)
#endif // ENABLE_SERIALIZATION
