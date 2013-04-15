// $Id$
#include <iostream>

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/mechanics/EngineeringStress3D.h"

// constructor
NuTo::EngineeringStress3D::EngineeringStress3D(): ConstitutiveOutputBase::ConstitutiveOutputBase()
{
    for (unsigned int count = 0; count < 6; count++)
    {
        (*this)[count] = 0.0;
    }
}

// number of components
unsigned int NuTo::EngineeringStress3D::GetNumberOfComponents() const
{
    return 6;
}

// get Engineering stress
const double* NuTo::EngineeringStress3D::GetData() const
{
    return data();
}

//! @brief ... sets the components of the Engineering stress tensor
//! @param ... components of the Engineering stress tensor (stored in an array)
//! @sa mEngineeringStress
void NuTo::EngineeringStress3D::SetData(double rData[6])
{
    for (unsigned int count = 0; count < 6; count++)
    {
    	(*this)[count] = rData[count];
    }
}


// info routine
void NuTo::EngineeringStress3D::Info(unsigned short rVerboseLevel) const
{
    std::cout << "    components of Engineering stress tensor (vector notation): "
              << (*this)[0] << ", " << (*this)[1] << ", " << (*this)[2] << ", "
              << (*this)[3] << ", " << (*this)[4] << ", " << (*this)[5] << std::endl;
}

#ifdef ENABLE_SERIALIZATION

//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::EngineeringStress3D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::EngineeringStress3D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::EngineeringStress3D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::EngineeringStress3D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::EngineeringStress3D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::EngineeringStress3D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::EngineeringStress3D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize EngineeringStress3D" << std::endl;
#endif
    ar &  BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveOutputBase);
    ar & boost::serialization::make_nvp ("EngineeringStrain3DEigen",boost::serialization::base_object< FullVectorFixed<double,6> > ( *this ));
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize EngineeringStress3D" << std::endl;
#endif
}
#endif // ENABLE_SERIALIZATION
