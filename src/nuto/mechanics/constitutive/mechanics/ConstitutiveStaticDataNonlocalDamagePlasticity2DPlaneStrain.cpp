// $ld: $ 
// ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain.cpp
// created May 6, 2010 by Joerg F. Unger


#include "nuto/mechanics/constitutive/mechanics/ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain.h"
#include "nuto/mechanics/constitutive/mechanics/ConstitutiveStaticDataPrevEngineeringStressStrain2DPlaneStrain.h"

//! @brief constructor
NuTo::ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain::ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain()
   : NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain2DPlaneStrain::ConstitutiveStaticDataPrevEngineeringStressStrain2DPlaneStrain()
{
	mKappa = 0.;

	mEpsilonP[0] = 0.;
	mEpsilonP[1] = 0.;
	mEpsilonP[2] = 0.;
	mEpsilonP[3] = 0.;
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template<class Archive>
void NuTo::ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain::serialize(Archive & ar, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveStaticDataPrevEngineeringStressStrain2DPlaneStrain)
       & BOOST_SERIALIZATION_NVP(mKappa)
       & BOOST_SERIALIZATION_NVP(mEpsilonP);
}
#endif // ENABLE_SERIALIZATION

//!@ brief reinterpret as nonlocal damage2d static data
NuTo::ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain* NuTo::ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain::AsNonlocalDamagePlasticity2DPlaneStrain()
{
    return this;
}

//!@ brief reinterpret as nonlocal damage2d static data
const NuTo::ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain* NuTo::ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain::AsNonlocalDamagePlasticity2DPlaneStrain()const
{
    return this;
}